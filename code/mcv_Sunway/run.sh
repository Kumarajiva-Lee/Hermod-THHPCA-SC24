bsub -J DSL -p -I -o DSL.log -q q_share -share_size 15000 -host_stack 1024 -cgsp 64 -b -n 6 -exclu ./app
