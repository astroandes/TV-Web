EXEC_PATH=/srv/cosmusr/scratch/jforero/TV-Web/src/
SNAP=088 #0062 #0085 #0052
NGRID=1024
DATA_PATH_IN=/srv/cosmdata/multidark/MD_3840_Planck1/SNAPS/snap_$SNAP
DATA_PATH_OUT=/srv/cosmusr/scratch/jforero/MDPL/T-Web/$NGRID/
mkdir -p $DATA_PATH_OUT



CIC_INTERPOLATION=1
TWEB_GAUSSIAN=0

if(($CIC_INTERPOLATION == 1)); then
    $EXEC_PATH/./smooth_bolshoi.x $DATA_PATH_IN/snap_$SNAP -n $NGRID -outdir $DATA_PATH_OUT 
fi



if(($TWEB_GAUSSIAN == 1));then
    smooth=(1.0)
    for i in "${smooth[@]}"
    do
        $EXEC_PATH./smooth_bolshoi.x -d $DATA_PATH_OUT/PMcrsFULL.$SNAP.DAT  -outdir $DATA_PATH_OUT -s $i
        $EXEC_PATH./smooth_bolshoi.x -d $DATA_PATH_OUT/PMcrsFULL.$SNAP.DAT  -outdir $DATA_PATH_OUT -eigenvector -s $i
    done
fi
