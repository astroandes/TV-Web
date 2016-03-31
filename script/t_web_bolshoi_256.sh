EXEC_PATH=/srv/cosmusr/scratch/jforero/TV-Web/src/
SNAP=0354
NGRID=256
DATA_PATH_IN=/srv/cosmdata/bolshoi/PMFILES/
DATA_PATH_OUT=/srv/cosmusr/scratch/jforero/Bolshoi/T-Web/$NGRID/
mkdir -p $DATA_PATH_OUT



CIC_INTERPOLATION=1
TWEB_GAUSSIAN=0


if(($CIC_INTERPOLATION == 1)); then
   snap_index_list=(0 1 2 3 4 5 6 7)
    for i in "${snap_index_list[@]}"
    do
    $EXEC_PATH/./smooth_bolshoi.x $DATA_PATH_IN/PMcrs$i.$SNAP.DAT -artfile $DATA_PATH_IN/PMcrd.$SNAP.DAT -artnfiles 16 -myart $i -n $NGRID -outdir $DATA_PATH_OUT
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
