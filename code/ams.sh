 
# for j in 1 2 3 4 5 6
# do
# nohup ./run -eta_0 0.08 > stdout.txt 2>stderr.txt &
# nohup ./run -eta_0 0.10 > stdout.txt 2>stderr.txt &
# nohup ./run -eta_0 0.12 > stdout.txt 2>stderr.txt &
# nohup ./run -eta_0 0.14 > stdout.txt 2>stderr.txt &
# nohup ./run -eta_0 0.16 > stdout.txt 2>stderr.txt &
# nohup ./run -eta_0 0.18 > stdout.txt 2>stderr.txt &
#local eta
# eta=$((0.08+$j*0.02))
./run -eta_0 0.12
./run -eta_0 0.12
./run -eta_0 0.12
./run -eta_0 0.12
# ./run -eta_0 0.08 -E_PCG 1 
# ./run -eta_0 0.10 -E_PCG 1 
# ./run -eta_0 0.12 -E_PCG 1
# ./run -eta_0 0.14 -E_PCG 1 
# ./run -eta_0 0.16 -E_PCG 1 
# ./run -eta_0 0.18 -E_PCG 1 
# ./run -eta_0 0.08 -E_PCG 0
# ./run -eta_0 0.10 -E_PCG 0 
# ./run -eta_0 0.12 -E_PCG 0 
# ./run -eta_0 0.14 -E_PCG 0 
# ./run -eta_0 0.16 -E_PCG 0 
# ./run -eta_0 0.18 -E_PCG 0 


# ./run -eta_0 0.08
./run -eta_0 0.10
./run -eta_0 0.12
#./home/cfeng/MCSS/a.out -eta_0 0.08
#./home/cfeng/MCSS/a.out -eta_0 0.10
#./home/cfeng/MCSS/a.out -eta_0 0.12
./run -eta_0 0.14
./run -eta_0 0.16
./run -eta_0 0.18

# done
