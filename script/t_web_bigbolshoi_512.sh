EXEC_PATH=/srv/cosmusr/scratch/jforero/TV-Web/src/
SNAP=0085 #0062 #0085 #0052
NGRID=512
DATA_PATH_IN=/srv/cosmdata/multidark/MDark_2048_om0.27/PMFILES/
DATA_PATH_OUT=/srv/cosmusr/scratch/jforero/MultiDark/T-Web/$NGRID/
mkdir -p $DATA_PATH_OUT



CIC_INTERPOLATION=0
TWEB_GAUSSIAN=1

if(($CIC_INTERPOLATION == 1)); then
    snap_index_list=(0 1 2 3 4 5 6 7 8 9)
    for i in "${snap_index_list[@]}"
    do
    $EXEC_PATH/./smooth_bolshoi.x $DATA_PATH_IN/PMcrs0$i.$SNAP.DAT -artfile $DATA_PATH_IN/PMcrd.$SNAP.DAT -artnfiles 16 -m\
yart $i -n $NGRID -outdir $DATA_PATH_OUT 
    done
    snap_index_list=(10 11 12 13 14 15)
    for i in "${snap_index_list[@]}"
    do
        $EXEC_PATH/./smooth_bolshoi.x $DATA_PATH_IN/PMcrs$i.$SNAP.DAT -artfile $DATA_PATH_IN/PMcrd.$SNAP.DAT -artnfiles 16\
 -myart $i -n $NGRID -outdir $DATA_PATH_OUT 
    done
fi



if(($TWEB_GAUSSIAN == 1));then
    smooth=(1.0)
    for i in "${smooth[@]}"
    do
        $EXEC_PATH./smooth_bolshoi.x -d $DATA_PATH_OUT/PMcrsFULL.$SNAP.DAT  -outdir $DATA_PATH_OUT -s $i
        $EXEC_PATH./smooth_bolshoi.x -d $DATA_PATH_OUT/PMcrsFULL.$SNAP.DAT  -outdir $DATA_PATH_OUT -eigenvector -s $i
    done
fi
