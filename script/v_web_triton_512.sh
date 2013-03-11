EXEC_PATH=/home/forero/Dropbox/TV-Web/src
NGRID=512
DATA_PATH_IN=/home/forero/work/data/Triton/snap036/snap
DATA_PATH_VWEB_OUT=/home/forero/work/data/Triton/snap036/vweb/$NGRID
DATA_PATH_TWEB_OUT=/home/forero/work/data/Triton/snap036/tweb/$NGRID
mkdir -p $DATA_PATH_VWEB_OUT
mkdir -p $DATA_PATH_TWEB_OUT
EXEC_NAME=smooth_bolshoi.x
SNAP_NAME=snapshot_036

CIC_INTERPOLATION=0 #makes the CIC interpolation                                                                                      
VWEB_CIC=0 # get the v-web for the CIC interpolation     
VWEB_GAUSSIAN=0 #get the v-web for the gaussian interpolation
TWEB_CIC=1 #old t-web for CIC interpolation
TWEB_GAUSSIAN=1 #old t-web for gaussian interpolation

if(($CIC_INTERPOLATION == 1)); then
    $EXEC_PATH/./$EXEC_NAME $DATA_PATH_IN/$SNAP_NAME -n $NGRID -outdir $DATA_PATH_VWEB_OUT -vel
fi

if(($VWEB_CIC == 1));then
    $EXEC_PATH/./$EXEC_NAME -d $DATA_PATH_VWEB_OUT/$SNAP_NAME -outdir $DATA_PATH_VWEB_OUT -vel 
    $EXEC_PATH/./$EXEC_NAME -d $DATA_PATH_VWEB_OUT/$SNAP_NAME  -outdir $DATA_PATH_VWEB_OUT -vel -eigenvector
fi

if(($VWEB_GAUSSIAN == 1));then
    smooth=(1.0 2.0)
    for i in "${smooth[@]}"
    do
	$EXEC_PATH/./$EXEC_NAME -d $DATA_PATH_VWEB_OUT/$SNAP_NAME -outdir $DATA_PATH_VWEB_OUT -vel -s $i
	$EXEC_PATH/./$EXEC_NAME -d $DATA_PATH_VWEB_OUT/$SNAP_NAME  -outdir $DATA_PATH_VWEB_OUT -vel -eigenvector -s $i
    done
fi

if(($TWEB_CIC == 1));then
        $EXEC_PATH/./$EXEC_NAME -d $DATA_PATH_TWEB_OUT/$SNAP_NAME  -outdir $DATA_PATH_TWEB_OUT 
        $EXEC_PATH/./$EXEC_NAME -d $DATA_PATH_TWEB_OUT/$SNAP_NAME  -outdir $DATA_PATH_TWEB_OUT -eigenvector 
fi

if(($TWEB_GAUSSIAN == 1));then
    smooth=(1.0 2.0)
    for i in "${smooth[@]}"
    do
        $EXEC_PATH/./$EXEC_NAME -d $DATA_PATH_TWEB_OUT/$SNAP_NAME  -outdir $DATA_PATH_TWEB_OUT -s $i
        $EXEC_PATH/./$EXEC_NAME -d $DATA_PATH_TWEB_OUT/$SNAP_NAME  -outdir $DATA_PATH_TWEB_OUT -eigenvector -s $i
    done
fi
