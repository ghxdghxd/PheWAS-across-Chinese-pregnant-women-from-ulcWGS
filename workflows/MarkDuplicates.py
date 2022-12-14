#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-07-27 12:55:09
# @Author  : JT Guo
# @Email   : guojt-4451@163.com
# @Version : 1.0

"""MarkDuplicates."""


import os
import re
import sys
import subprocess
import tempfile
import multiprocessing
# import glob

version = "3.0"


def get_args():
    """Get arguments from commond line"""
    try:
        import argparse
    except ImportError as imerr:
        print("\033[1;31m" + str(imerr) + " \033[0m")
        sys.exit()

    parser = argparse.ArgumentParser(
        usage="%(prog)s", fromfile_prefix_chars='@', description=__doc__)

    parser.add_argument("--path2GATK",
                        metavar="PATH",
                        dest="GATK",
                        required=True,
                        help="path to gatk")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i",
                       metavar="FILE",
                       dest="bam_file",
                       help="input bamfile")
    group.add_argument("-I",
                       metavar="FILE",
                       dest="bam_list",
                       help="list of input bamfile")

    parser.add_argument("-o",
                        metavar="DIR",
                        dest="output_dir",
                        default=os.environ['HOME'],
                        help="output dir or output file \
                        [" + os.environ['HOME'] + "]")
    parser.add_argument("-p",
                        metavar="INT",
                        default=1,
                        dest="processes_number",
                        type=int,
                        help="analyze multiple samples simultaneously [1]")
    parser.add_argument("-t",
                        metavar="INT",
                        default=1,
                        dest="threads_number",
                        type=int,
                        help="number of threads to allocate \
                        to each sample [1]")
    parser.add_argument("-m",
                        metavar="STR",
                        default="2g",
                        dest="Xmx",
                        help="memory [2g]")

    parser.add_argument("--qsub",
                        action="store_true",
                        default=False,
                        dest="qsub",
                        help="run crest in cluster")
    parser.add_argument("--remove_duplicates",
                        action="store_true",
                        default="false",
                        dest="rmdup",
                        help="remove PCR duplicates")
    parser.add_argument("--validation_stringency",
                        type=str,
                        default="STRICT",
                        metavar="{STRICT, LENIENT, SILENT}",
                        dest="validation_stringency",
                        help="Validation stringency for all SAM files \
                        read by this program [STRICT]")
    parser.add_argument("--node",
                        metavar="STR",
                        dest="node",
                        help="name of nodes")
    parser.add_argument("-n",
                        metavar="INT",
                        default=1,
                        dest="nodes_number",
                        type=int,
                        help="number of nodes [1]")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    else:
        args = parser.parse_args()
        args.output_dir = os.path.realpath(args.output_dir)
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

def mark_duplicates(bam, run = True):
    bam_name = os.path.basename(bam).split(".sorted.bam")[0]
    cmd = args.GATK + " MarkDuplicates"
    cmd += " --REMOVE_DUPLICATES " + str(args.rmdup)
    cmd += " --VALIDATION_STRINGENCY " + args.validation_stringency
    cmd += " --INPUT " + bam
    cmd += " --OUTPUT " + args.output_dir + "/" + bam_name + ".sorted.dedup.bam"
    cmd += " --METRICS_FILE " + args.output_dir + "/" + bam_name + ".metrics"
    cmd += " --CREATE_INDEX true --ASSUME_SORTED true"
    if run:
        print(cmd)
        os.system(cmd)
    else:
        return bam_name, cmd

def qsub_mark_duplicates(bam):
    bam_name, cmd = mark_duplicates(bam, run=False)
    ftmp = tempfile.NamedTemporaryFile()
    ftmp.write(b"#!/bin/bash\n")
    ftmp.write(b"#PBS -N markdup-" + bam_name.encode() + b"\n")
    ftmp.write(b"#PBS -o " + args.output_dir.encode() +
               b"/" + bam_name.encode() + b".markdup\n")
    if args.node:
        ftmp.write(b"#PBS -l nodes=1:" + args.node.encode() +
                   b":ppn=" + args.threads_number.encode() +
                   b",mem=" + args.Xmx.encode() + b",walltime=100:00:00\n")
    else:
        ftmp.write(b"#PBS -l nodes=1:ppn=" + args.threads_number.encode() +
                   b",mem=" + args.Xmx.encode() + b",walltime=100:00:00\n")
    ftmp.write(b"#PBS -j oe\ncd $PBS_O_WORKDIR\nsource /etc/profile.d/set.sh\n")
    ftmp.write(cmd)
    ftmp.seek(0)
    print(ftmp.read())
    # os.system("qsub " + ftmp.name)
    ftmp.close()

def main():
    global args
    args = get_args()
    pool = multiprocessing.Pool(processes=int(args.processes_number))
    
    if args.bam_list:
        with open(args.bam_list) as f:
            bam_list = map(lambda x: x.strip(), f.readlines())
            if args.qsub:
                pool.map(qsub_mark_duplicates, bam_list)
            else:
                pool.map(mark_duplicates, bam_list)
            pool.close()
            pool.join()
    elif args.bam_file:
        if args.qsub:
            qsub_mark_duplicates(args.bam_file)
        else:
            mark_duplicates(args.bam_file)

if __name__ == '__main__':
    sys.exit(main())


