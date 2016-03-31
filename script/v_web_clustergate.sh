EXEC_PATH=/lustre/home/ciencias/fisica/je.forero/TV-Web/src/
NGRID=32
DATA_PATH_IN=/lustre/home/ciencias/fisica/je.forero/Gadget-2.0.7/Gadget2/lcdm_ancos/
DATA_PATH_VWEB_OUT=/lustre/home/ciencias/fisica/je.forero/test/$NGRID
DATA_PATH_TWEB_OUT=/lustre/home/ciencias/fisica/je.forero/test/$NGRID
mkdir -p $DATA_PATH_VWEB_OUT
mkdir -p $DATA_PATH_TWEB_OUT
#/lustre/home/ciencias/fisica/je.forero/BOGOTA_School/Box256/
EXEC_NAME=smooth_bolshoi.x
SNAP_NAME=snapshot_005

CIC_INTERPOLATION=1 #makes the CIC interpolation                                                                                      
VWEB_CIC=1 # get the v-web for the CIC interpolation     
VWEB_GAUSSIAN=1 #get the v-web for the gaussian interpolation
TWEB_CIC=0 #old t-web for CIC interpolation
TWEB_GAUSSIAN=0 #old t-web for gaussian interpolation

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
