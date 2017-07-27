#!/bin/sh

#===========================================================================
#*        AUTO-SUBMISSION THE JOBS FOR CALCULATING I-V CURVE               *
#*                                                                         *  
#*                        -- BY CHEN SHUGUANG                              * 
#*                          APRIL 30TH, 2015                               * 
#===========================================================================


#----------------------------       parameters       -----------------------

ngrid=18;           # number of voltage points need to be calculated
bot=-0.8;           # The lowest bias voltage (Volt = right - left)
top=1.0;            # The highest bias voltage
nnode=6;            # The smallest number of OMP cores requested

echo "Total voltage points number: " $ngrid;
echo "lowest bias voltage: " $bot;
echo "Highest bias voltage: " $top;
echo "OMP nodes needed: "$nnode;

#===========================================================================

for i in $(seq 0 $ngrid)

do
  echo "$[i+1]-th job";
  volt=$(echo |awk '{print '"(${top}"' - '"(${bot}))"' * '"$i"' / '"$ngrid"' + '"($bot)"'}');
  # Calculate the voltage in the i-th job ( Note: Variales in AWK should be cited as '"$VAR"' )

  lvolt=$(echo |awk '{print  '"-($volt)"'/2}');
  rvolt=$(echo |awk '{print  '"$volt"'/2}');
  # Now the bias is applied symmytrically 


  echo "Voltage applied : " $volt;
  echo "left side bias : " $lvolt;
  echo "right side bias: " $rvolt;


  dname=b${volt}.gate;
  # The name of the directory for i-th job

  mkdir $dname;
  cp -r skf $dname;
  cp in.gstd qsub.td lead*.dat gridpoint.xyz metal.in  $dname;
  # Copy the input files needed to the working directory

  cd ${dname};
  # Enter the working directory 

  nline=$(awk '/\$transport/ {print NR}' in.gstd);
  # Find the line in the input file containning the parameters of voltage

  sed -i "${nline}s/voltage(1)=-0.05d0/voltage(1)=${lvolt}/" in.gstd;
  sed -i "${nline}s/voltage(2)=0.05d0/voltage(2)=${rvolt}/"  in.gstd;
  # replace the original voltage (now -0.05d0/0.05d0) with the one we need

  while true;
  do

    qstat -f |awk -F'[ /]+' '/all.q@compute-1/ && $8 == ""  {print $1,$2,$3,$4,$5,$6,$7,$8}' \
    |sort -n -k 4\
    |awk -F'[ ./]+' 'NR==1 {print $1 "." $2, $7-$6}' >status.temp;
    # List the occupation status of the clusters; 
    # List those nodes belongs to GHC(all.q@C1-*) and are working ($8 may be "d" or "AU" when dead);
    # Sort the nodes available by their occupations and select the one with the most free cores; 
    # Store the node name as well as the number of the cores available in a temporary file
    
    inode=$(awk '{print $1}' status.temp);
    ncore=$(awk '{print $2}' status.temp);
    # Extract the node name
    # Extract the number of the free cores
    
    rm status.temp;
    # Delete the temporary file 
    
    echo "The nodes with lowest loads availble: " $inode;
    echo "Free core availble: " $ncore;

    if [ "$ncore" -ge "$nnode" ]; then
      break;
    else
      sleep 30s;
    fi
    # Check whether there is any machine satisfy our needs (having enough cores)
    # If yes, continue to submit the job; Elsewise wait for 30s and startover 

  done  

  nline=$(awk '/pe orte/ {print NR}' qsub.td);
  sed -i "${nline}s/pe orte .*/pe orte ${ncore}/" qsub.td;
  # Find the lines in the qsub script containning the OMP settings
  # Assign the core number to be the largest one we have


  nline=$(awk '/H-ppc.nogate/ {print NR}' qsub.td);
  sed -i "${nline}s/H-ppc.nogate/b${volt}.gate/" qsub.td;
  # Find the line in the qsub script containing the mission name
  # Assign a new name containning the bias

  qsub -q "$inode" qsub.td;
  # Submit the mission

  cd ..;
  # Exit the working directory (preparing for the next job)

  sleep 20s;
  # Wati 20s before the next job starts, avoiding conflict

done


