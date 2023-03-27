from time import sleep
from celery import shared_task
from django.core.mail import send_mail, mail_admins, BadHeaderError
import mstid.music_support as msc
import datetime
import os

@shared_task
def notify_customers(current_user_email,current_user_isauth):
    print("Sending 10k emails...")
    # print(message)
    sleep(1)
    print("Emails were successfully sent!")
    if current_user_isauth:
        try:
            send_mail("Notification for completion of MSTID Classification","Classification has be completed","fransway26@outlook.com",[current_user_email])
        except BadHeaderError:
            pass


@shared_task
def run_classification(dct, mstid_index, new_list, 
                       recompute, reupdate_db, music_process, 
                       music_new_list, music_reupdate_db,
                       nprocs, multiproc,current_user_email,current_user_isauth):
    import matplotlib
    matplotlib.use('Agg')
    import mstid
    from mstid import run_helper
    
    # List start date
    list_start_year = int(dct["list_sDate"][:4])
    list_start_month = int(dct["list_sDate"][5:7])
    list_start_day = int(dct["list_sDate"][8:10])
    # List end date
    list_end_year = int(dct["list_eDate"][:4])
    list_end_month = int(dct["list_eDate"][5:7])
    list_end_day = int(dct["list_eDate"][8:10])
    
    dct["list_sDate"] = datetime.datetime(list_start_year,list_start_month,list_start_day)
    dct["list_eDate"] = datetime.datetime(list_end_year,list_end_month,list_end_day)
    db_name          = 'mstid'
    dct_list                    = run_helper.create_music_run_list(**dct)
    # Classification parameters go here. ###########################################
    classification_path = 'mstid_data/classification'

    #******************************************************************************#
    # No User Input Below This Line ***********************************************#
    #******************************************************************************#

    if mstid_index:
        # import ipdb;ipdb.set_trace()
        # Generate MSTID List and do rti_interp level processing.
        run_helper.get_events_and_run(dct_list,process_level='rti_interp',new_list=new_list,
                recompute=recompute,multiproc=multiproc,nprocs=nprocs)

        # Reload RTI Data into MongoDb. ################################################
        if reupdate_db:
            for dct in dct_list:
                mstid.updateDb_mstid_list(multiproc=multiproc,nprocs=nprocs,**dct)

        for dct in dct_list:
            # Determine if each event is good or bad based on:
            #   1. Whether or not data is available.
            #   2. Results of pyDARNmusic.utils.checkDataQuality()
            #       (Ensures radar is operational for a minimum amount of time during the data window.
            #        Default is to require the radar to be turned off no more than 10 minutes in the
            #        data window.)
            #   3. The fraction of radar scatter points present in the data window.
            #       (Default requires minimum 67.5% data coverage.)
            #   4. The percentage of daylight in the data window.
            #       (Default requires 100% daylight in the data window.)
            mstid.classify.classify_none_events(**dct) 

            # Generate a web page and copy select figures into new directory to make it easier
            # to evaluate data and see if classification algorithm is working.
            mstid.classify.rcgb(classification_path=classification_path,**dct)

        # Run FFT Level processing on unclassified events.
        run_helper.get_events_and_run(dct_list,process_level='fft',category='unclassified',
                multiproc=multiproc,nprocs=nprocs)

        # Now run the real MSTID classification.
        mstid.classify.run_mstid_classification(dct_list,classification_path=classification_path,
                multiproc=multiproc,nprocs=5)

        print('Plotting calendar plot...')
        mstid.calendar_plot(dct_list,db_name=db_name)


    # Run actual MUSIC Processing ##################################################
    
    if music_process:
        for dct in dct_list:
            dct['input_mstid_list']     = dct['mstid_list']
            dct['input_db_name']        = dct['db_name']
            dct['input_mongo_port']     = dct['mongo_port']
            dct['mstid_list']           = 'music_'+dct['mstid_list']
            dct['data_path']            = 'mstid_data/music_data'
            dct['hanning_window_space'] = True
    #        dct['bad_range_km']         = 500 # Set to 500 for MUSIC Calculation
            dct['bad_range_km']         = None # Set to None to match original calculations
        
        run_helper.get_events_and_run(dct_list,process_level='music',
                new_list=music_new_list,category=['mstid','quiet'],
                multiproc=multiproc,nprocs=nprocs)

        if music_reupdate_db:
            for dct in dct_list:
                mstid.updateDb_mstid_list(multiproc=multiproc,nprocs=nprocs,**dct)

    #tunnel.kill()
    print("I'm done!")
    if current_user_isauth:
        try:
            send_mail("Notification for completion of MSTID Classification","Classification has be completed","fransway26@outlook.com",[current_user_email])
        except BadHeaderError:
            pass

@shared_task
def create_music_object_task(radar, sDatetime, fDatetime,
                             beamLimits_0, beamLimits_1, gateLimits_0,
                             gateLimits_1, interpolationResolution,
                             filterNumtaps, firFilterLimits_0, 
                             firFilterLimits_1, window_data,
                             kx_max,ky_max, autodetect_threshold_str,
                             neighborhood_0, neighborhood_1,
                             current_user_email, current_user_isauth):
    
    try:
        bl0 = int(beamLimits_0)
    except:
        bl0 = None
    try:
        bl1 = int(beamLimits_1)
    except:
        bl1 = None
    beamLimits = (bl0, bl1)

    try:
        gl0 = int(gateLimits_0)
    except:
        gl0 = None
    try:
        gl1 = int(gateLimits_1)
    except:
        gl1 = None
    gateLimits = (gl0,gl1)

    try:
        interpRes = int(interpolationResolution)
    except:
        interpRes = None

    try:
        numtaps = int(filterNumtaps)
    except:
        numtaps = None
    
    try:
        cutoff_low  = float(firFilterLimits_0)
    except:
        cutoff_low  = None

    try:
        cutoff_high  = float(firFilterLimits_1)
    except:
        cutoff_high  = None

    try:
        kx_max  = float(kx_max)
    except:
        kx_max  = 0.05

    try:
        ky_max  = float(ky_max)
    except:
        ky_max  = 0.05

    try:
        autodetect_threshold  = float(autodetect_threshold_str)
    except:
        autodetect_threshold  = 0.35
    
    try:
        nn0 = int(neighborhood_0)
    except:
        nn0 = None
    try:
        nn1 = int(neighborhood_1)
    except:
        nn1 = None
    neighborhood = (nn0,nn1)

    ################################################################################ 

    musicPath   = msc.get_output_path(radar, sDatetime, fDatetime)

    try:
        import shutil
        shutil.rmtree(musicPath)
    except:
        pass

    if window_data == 'true':
        window_data = True
    else:
        window_data = False
        

    dataObj = msc.createMusicObj(radar.lower(), sDatetime, fDatetime
        ,beamLimits                 = beamLimits
        ,gateLimits                 = gateLimits
        ,interpolationResolution    = interpRes
        ,filterNumtaps              = numtaps 
        ,fitfilter                  = True
        )
    # import ipdb;ipdb.set_trace()

    picklePath  = msc.get_pickle_name(radar,sDatetime,fDatetime,getPath=True,createPath=False)


    # Create a run file. ###########################################################
    runParams = {}
    runParams['radar']              = radar.lower()
    runParams['sDatetime']          = sDatetime
    runParams['fDatetime']          = fDatetime
    runParams['beamLimits']         = beamLimits
    runParams['gateLimits']         = gateLimits
    runParams['interpRes']          = interpRes
    runParams['filter_numtaps']     = numtaps
    runParams['filter_cutoff_low']  = cutoff_low
    runParams['filter_cutoff_high'] = cutoff_high
    runParams['path']               = musicPath
    runParams['musicObj_path']      = picklePath
    runParams['window_data']        = window_data
    runParams['kx_max']             = kx_max
    runParams['ky_max']             = ky_max
    runParams['autodetect_threshold'] = autodetect_threshold
    runParams['neighborhood']        = neighborhood
  
    msc.Runfile(radar.lower(), sDatetime, fDatetime, runParams)

    # Generate general RTI plot for original data. #################################
    #Mark where sampling window starts and stops.
    dataObj.DS000_originalFit.metadata['timeLimits'] = [runParams['sDatetime'],runParams['fDatetime']]

    rti_beams   = msc.get_default_beams(runParams,dataObj)
    rtiPath     = os.path.join(musicPath,'000_originalFit_RTI.png')

    msc.plot_music_rti(dataObj,fileName=rtiPath,dataSet="originalFit",beam=rti_beams)
     
    if current_user_isauth:
        try:
            send_mail("Notification for completion of MUSIC Object Creation","MUSIC Object has been created","fransway26@outlook.com",[current_user_email])
        except BadHeaderError:
            pass  

@shared_task
def run_music_task(runfile_path,current_user_email, current_user_isauth):
    msc.run_music(runfile_path)
    msc.music_plot_all(runfile_path)
    
    if current_user_isauth:
        try:
            send_mail("Notification for completion of MSTID Analysis","MSTID Analysis has be completed","fransway26@outlook.com",[current_user_email])
        except BadHeaderError:
            pass