#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Version : 1.0

"""Base_quality_score_recalibration."""


import os
import re
import sys
import subprocess
import tempfile
# import glob

version = "3.0"


def get_args3():
    """Get arguments from commond line."""
    try:
        import argparse
    except ImportError as imerr:
        print("\033[1;31m" + str(imerr) + " \033[0m")
        sys.exit()

    parser = argparse.ArgumentParser(prog="runGATK4",
                                     usage="%(prog)s",
                                     fromfile_prefix_chars='@',
                                     description=__doc__)

    parser.add_argument('-v', '--version',
                        action='version',
                        version="%(prog)s " + version)

    parser.add_argument("--path2GATK",
                        metavar="PATH",
                        dest="GATK",
                        required=True,
                        help="path to gatk")

    parser.add_argument("-r",
                        metavar="FILE",
                        dest="ref_fa",
                        required=True,
                        help="faidx indexed reference sequence file")
    parser.add_argument("-i",
                        metavar="FILE",
                        dest="input_list",
                        action='append',
                        required=True,
                        help="list of input lane_bam_files \
                            (one lane_bam per line)")
    parser.add_argument("-o",
                        metavar="FILE",
                        dest="output_file",
                        help="output file name without extension: .bam")
    parser.add_argument("-l",
                        metavar="FILE",
                        dest="intervals",
                        help="One or more genomic intervals over which to operate")
    parser.add_argument("--known-sites",
                        metavar="FILE",
                        action='append',
                        dest="knownsites",
                        help="knownsites")
    parser.add_argument("--java-Xmx",
                        metavar="STR",
                        dest="java_Xmx",
                        default="4G",
                        help="memory [4g]")

    spark_group = parser.add_argument_group("spark arguments")
    spark_group.add_argument("--spark",
                                action="store_true",
                                default=False,
                                dest="spark",
                                help="run in spark pipeline [False]")
    spark_group.add_argument("-t",
                            metavar="INT",
                            default=1,
                            dest="threads_number",
                            type=int,
                            help="number of threads/per_sample using Spark[1]")

    clustered_group = parser.add_argument_group("Clustered arguments")
    clustered_group.add_argument("--qsub",
                                 action="store_true",
                                 default=False,
                                 dest="qsub",
                                 help="run crest in cluster [False]")
    clustered_group.add_argument("--node",
                                 metavar="STR",
                                 dest="node",
                                 type=str,
                                 help="name of nodes (e.g: n1,n2,...)")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        # check_dependencies(["group_name"])
        args = parser.parse_args()
        args.threads_number = str(args.threads_number)
        return args


def check_dependencies(tools):
    """Ensure required tools are present."""
    if not isinstance(tools, list):
        tools = [tools]
    try:
        for tool in tools:
            try:
                subprocess.check_output(["which", tool]).strip()
            except subprocess.CalledProcessError:
                print("\033[1;31m" + __file__ + " requires " +
                      tool + "\033[0m")
    except:
        sys.exit()


def makedir(new_dir, exist_dir=None):
    """Make a directory.

    If it doesn't exist, handling concurrent race conditions.
    """
    if exist_dir:
        new_dir = os.path.join(exist_dir, new_dir)
    if os.path.exists(new_dir):
        print("The " + new_dir + " is already exist")
    else:
        print("Make " + new_dir)
        os.makedirs(new_dir)


def get_same_name(str1, str2):
    """Return the same pattern between str1 and str2."""
    for i in range(1, len(str1)):
        if str1[:-i] == str2[:-i] and str1[:-i][-1].isalnum():
            return os.path.basename(str1[:-i])


# def merge_lane(args):
#     """Picard Merge Sam Files"""
#     merge_cmd = args.GATK + " MergeSamFiles" 
#     merge_cmd += " -I " + " -I ".join(args.input_list)
#     merge_cmd += " -O " + args.output_file + ".merge.bam"
#     merge_cmd += " --java-options \"-Xmx" + args.java_Xmx + \
#         " -Djava.io.tmpdir=" + os.path.dirname(args.output_file) + "/tmp\""
#     return merge_cmd

def bqsr_pipeline_spark(bam, args):
    """BQSRPipelineSpark"""
    makedir(os.path.dirname(args.output_file) + "/tmp")
    cmd = args.GATK + " BQSRPipelineSpark"
    cmd += " -R " + args.ref_fa
    cmd += " -I " + bam
    cmd += " -O " + args.output_file + ".table"
    if os.path.exists(args.intervals):
        cmd += " -L " + args.intervals
    cmd += " --known-sites " + " --known-sites ".join(args.knownsites)
    cmd += " --spark-master local[" + args.threads_number + "]"
    cmd += " --java-options \"-Xmx" + args.java_Xmx + "\""
    cmd += " --tmp-dir " + os.path.dirname(args.output_file) + "/tmp"
    return cmd

def bqsr(bam, args):
    """Base_quality_score_recalibration."""
    makedir(os.path.dirname(bam))
    cmd = args.GATK + " BaseRecalibrator"
    cmd += " -R " + args.ref_fa
    cmd += " -I " + bam
    cmd += " -O " + args.output_file + ".table"
    if os.path.exists(args.intervals):
        cmd += " -L " + args.intervals
    cmd += " --known-sites " + " --known-sites ".join(args.knownsites)
    # cmd += " --spark-master local[" + args.threads_number + "]"
    cmd += " --java-options \"-Xmx" + args.java_Xmx + "\""
    cmd += " --tmp-dir " + os.path.dirname(args.output_file) + "/tmp"
    return cmd


def apply_bqsr(bam, args):
    """ApplyBQSR"""
    cmd = args.GATK + " ApplyBQSR"
    cmd += " -R " + args.ref_fa
    cmd += " -I " + bam
    cmd += " --bqsr-recal-file " + args.output_file + ".table"
    cmd += " -O " + args.output_file + ".bqsr.bam"
    if os.path.exists(args.intervals):
        cmd += " -L " + args.intervals
    cmd += " --java-options \"-Xmx" + args.java_Xmx + "\""
    cmd += " --tmp-dir " + os.path.dirname(args.output_file) + "/tmp"
    return cmd


def qsub(cmd_list, args, job_name):
    """Qsub."""
    with tempfile.NamedTemporaryFile() as ftmp:
        ftmp.write(b"#!/bin/bash\n")
        ftmp.write(b"#PBS -N " + job_name.encode() + b"\n")
        ftmp.write(b"#PBS -o " + args.output_file.encode() + b".log\n")
        if args.node:
            ftmp.write(b"#PBS -l nodes=1:" + args.node.encode() +
                    b":ppn=" + args.threads_number.encode() + b",mem=" + args.java_Xmx.encode() + b"\n")
        else:
            ftmp.write(b"#PBS -l nodes=1:ppn=" +
                    args.threads_number.encode() + b",mem=" + args.java_Xmx.encode() + b"\n")
        ftmp.write(b"#PBS -j oe\ncd $PBS_O_WORKDIR\n")
        ftmp.write(b"source /etc/profile.d/set.sh\n")
        # ftmp.write(b"rm " + os.path.dirname(args.output_file).encode() + b"/tmp/*" + b"\n")
        ftmp.write(cmd_list.encode())
        ftmp.seek(0)
        print(ftmp.read())
        # os.system("qsub " + ftmp.name)


def main():
    """Main."""
    args = get_args3()
    
    # if len(args.input_list) > 1:
    #     cmd_list = "\n".join([merge_lane(args), bqsr(
    #         args.output_file + ".merge.bam", args), apply_bqsr(
    #         args.output_file + ".merge.bam", args)])
    # else:
    #     # cmd_list = apply_bqsr(args.input_list[0], args)
    args.input_list[0] = os.path.realpath(args.input_list[0])
    
    if args.spark:
        cmd_list = bqsr_pipeline_spark(args.input_list[0], args)
    else:
        cmd_list = "\n".join([bqsr(args.input_list[0], args),
                        apply_bqsr(args.input_list[0], args)])

    # print(cmd_list)
    if args.qsub:
        qsub(cmd_list, args, "BQSR")
    else:
        os.system(cmd_list)

if __name__ == '__main__':
    sys.exit(main())


