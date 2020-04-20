#!/usr/bin/env bash
cloc --report-file=cloc_report.txt --list-file=./include.txt --exclude-list-file=./exclude.txt
