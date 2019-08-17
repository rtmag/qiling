python /root/qiling/juicer/misc/generate_site_positions.py \
Arima mm10_noScaffold /root/qiling/juicer/references/mm10_noScaffold.fasta

scp -r /root/qiling/juicer/restriction_sites e0056363@172.25.192.196:/hpctmp/e0056363/juicer

# In GPU cluster ssh e0056363@172.25.192.196
echo "bash /hpctmp/e0056363/juicer/scripts/juicer.sh \
-y /root/qiling/juicer/restriction_sites/mm10_noScaffold_Arima.txt \
-z /root/qiling/juicer/references/mm10_noScaffold.fasta -p chromsize"

qsub qiling_juicer.bash
