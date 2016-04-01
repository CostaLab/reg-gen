from __future__ import print_function
from __future__ import division
import urllib2
import sys

ALIAS_hg19="../data/hg19/alias.txt"
ALIAS_hg38="../data/hg38/alias.txt"
ALIAS_mm9="../data/mm9/alias.txt"

# HG19 and HG38
response = urllib2.urlopen('http://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=gd_prev_sym&col=gd_aliases&col=gd_pub_eg_id&col=gd_pub_refseq_ids&col=md_ensembl_id&status=Approved&status=Entry+Withdrawn&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit')
table = response.read()
# print(type(table))
alias = []

for line in table.split("\n"):
    l = line.split("\t")
    if l[-1]:
        # print(l)
        extra = "&".join([x for x in l if x])
        extra = extra.replace(", ", "&")
        # print(l[-1])
        alias.append("\t".join([l[-1],l[1],extra]))
        # sys.exit(0)
# print(len(table))

with open(ALIAS_hg19,"w") as f:
    for line in alias:
        print(line, file=f)

# with open(ALIAS_hg38,"w") as f:
#     for line in alias:
#         print(line, file=f)