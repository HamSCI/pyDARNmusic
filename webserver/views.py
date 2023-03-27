from django.shortcuts import render, redirect
from django.http import JsonResponse
from django.contrib import messages
from matplotlib import pyplot as plt 
from bson.objectid import ObjectId
from io import BytesIO
import matplotlib
matplotlib.use('Agg')
import datetime
import base64
import os,flash,jsonify
import pymongo
import pickle
import numpy as np
import glob
from matplotlib.figure import Figure

from pyDARNmusic import load_fitacf,music,musicRTP,checkDataQuality, stringify_signal_list,add_signal,del_signal
from mstid import mongo_tools
from mstid import more_music as mm
from mstid.stats_support import *
# import manual_support as mans
import mstid.music_support as msc
# Create your views here.

mongo         = pymongo.MongoClient()
db            = mongo.mstid

def send_email(request):
    from .tasks import notify_customers
    current_user = request.user
    current_user_email = "francistholley360@gmail.com"
    current_user_isauth = current_user.is_authenticated
    
    notify_customers(current_user_email,current_user_isauth)
    return render(request,"apphome.html")
#App HomePage
def apphome(request):
    return render(request,'apphome.html')

def home(request):
    return render(request,"home.html")

def classify_mstids(request):
    try:
        if request.POST:
            
            radars = []

            if request.POST.get("ade") != None :
                radars.append(request.POST.get("ade"))
            if request.POST.get("adw") != None : 
                radars.append(request.POST.get("adw"))
            if request.POST.get("bks") != None : 
                radars.append(request.POST.get("bks"))
            if request.POST.get("cve") != None :        
                radars.append(request.POST.get("cve"))
            if request.POST.get("cvw") != None : 
                radars.append(request.POST.get("cvw"))
            if request.POST.get("cly") != None : 
                radars.append(request.POST.get("cly"))
            if request.POST.get("fhe") != None : 
                radars.append(request.POST.get("fhe"))
            if request.POST.get("fhw") != None : 
                radars.append(request.POST.get("fhw"))
            if request.POST.get("gbr") != None : 
                radars.append(request.POST.get("gbr"))
            if request.POST.get("han") != None : 
                radars.append(request.POST.get("han"))
            if request.POST.get("hok") != None : 
                radars.append(request.POST.get("hok"))
            if request.POST.get("hkw") != None : 
                radars.append(request.POST.get("hkw"))
            if request.POST.get("inv") != None : 
                radars.append(request.POST.get("inv"))
            if request.POST.get("jme") != None : 
                radars.append(request.POST.get("jme"))
            if request.POST.get("kap") != None : 
                radars.append(request.POST.get("kap"))
            if request.POST.get("ksr") != None : 
                radars.append(request.POST.get("ksr"))
            if request.POST.get("kod") != None : 
                radars.append(request.POST.get("kod"))
            if request.POST.get("lyr") != None : 
                radars.append(request.POST.get("lyr"))
            if request.POST.get("pyk") != None : 
                radars.append(request.POST.get("pyk"))
            if request.POST.get("pgr") != None : 
                radars.append(request.POST.get("pgr"))
            if request.POST.get("rkn") != None : 
                radars.append(request.POST.get("rkn"))
            if request.POST.get("sas") != None : 
                radars.append(request.POST.get("sas"))
            if request.POST.get("sch") != None : 
                radars.append(request.POST.get("sch"))
            if request.POST.get("sto") != None : 
                radars.append(request.POST.get("sto"))
            if request.POST.get("wal") != None : 
                radars.append(request.POST.get("wal"))
            # Radars in Southern Hemisphere
            if request.POST.get("bpk") != None : 
                radars.append(request.POST.get("bpk"))
            if request.POST.get("dce") != None : 
                radars.append(request.POST.get("dce"))
            if request.POST.get("dcn") != None : 
                radars.append(request.POST.get("dcn"))
            if request.POST.get("fir") != None : 
                radars.append(request.POST.get("fir"))
            if request.POST.get("hal") != None : 
                radars.append(request.POST.get("hal"))
            if request.POST.get("ker") != None : 
                radars.append(request.POST.get("ker"))
            if request.POST.get("mcm") != None : 
                radars.append(request.POST.get("mcm"))
            if request.POST.get("san") != None : 
                radars.append(request.POST.get("san"))
            if request.POST.get("sps") != None : 
                radars.append(request.POST.get("sps"))
            if request.POST.get("sye") != None : 
                radars.append(request.POST.get("sye"))
            if request.POST.get("sys") != None : 
                radars.append(request.POST.get("sys"))
            if request.POST.get("tig") != None : 
                radars.append(request.POST.get("tig"))
            if request.POST.get("unw") != None : 
                radars.append(request.POST.get("unw"))
            if request.POST.get("zho") != None : 
                radars.append(request.POST.get("zho"))
            
            if len(radars) == 0:
                messages.add_message(request, messages.ERROR, 'At least one radar must be selected for classification')
                return render(request,"classify.html") 
            
            db_name          = 'mstid'
            dct                         = {}
            dct['radars']               = radars
            dct['db_name']              = db_name
            dct['data_path']            = 'mstid_data/mstid_index'
            
            # List start date
            list_start_year = int(request.POST.get("start_date")[:4])
            list_start_month = int(request.POST.get("start_date")[5:7])
            list_start_day = int(request.POST.get("start_date")[8:10])
            # List end date
            list_end_year = int(request.POST.get("end_date")[:4])
            list_end_month = int(request.POST.get("end_date")[5:7])
            list_end_day = int(request.POST.get("end_date")[8:10])
            if list_start_year >= list_end_year and list_start_month >= list_end_month and list_start_day >= list_end_day:
                messages.add_message(request, messages.ERROR, 'List Start Date must be smaller than List End Date')
                return render(request,"classify.html")
            
            dct['list_sDate']           = datetime.datetime(list_start_year,list_start_month,list_start_day)
            dct['list_eDate']           = datetime.datetime(list_end_year,list_end_month,list_end_day)

            dct['hanning_window_space'] = False
            if request.POST.get("hanning_window_space") != None:
                dct['hanning_window_space'] = True
            

            dct['bad_range_km'] = None
            if request.POST.get("bad_range_km"):
                if int(request.POST.get("bad_range_km")) < 0:
                    messages.add_message(request, messages.ERROR, 'bad_range_km must be a positive integer')
                    return render(request,"classify.html")
                dct['bad_range_km'] = request.POST.get("bad_range_km")
            
            
            mstid_index         = False
            if request.POST.get("mstid_index") != None:
                mstid_index = True
                
            if request.POST.get("mstid_index") != None and request.POST.get("hanning_window_space") != None:
                messages.add_message(request, messages.ERROR, 'Hanning Window Space must be unselected when Calculate MSTID Index is selected')
                return render(request,"classify.html")
            
            if request.POST.get("bad_range_km") and request.POST.get("mstid_index") != None:
                messages.add_message(request, messages.ERROR, 'bad_range_km must not have a value when Calculate MSTID Index is selected')
                return render(request,"classify.html")
            
            new_list = False
            if request.POST.get("new_list") != None:
                new_list = True

            recompute = False
            if request.POST.get("recompute") != None:
                recompute = True

            reupdate_db = False
            if request.POST.get("reupdate_db") != None:
                reupdate_db = True
            
            music_process = False
            if request.POST.get("music_process") != None:
                music_process = True
            
            music_new_list = False
            if request.POST.get("music_new_list") != None:
                music_new_list = True
            
            music_reupdate_db = False
            if request.POST.get("music_reupdate_db") != None:
                music_reupdate_db = True
            
            nprocs = 8
            if request.POST.get("nprocs"):
                nprocs = int(request.POST.get("nprocs"))

            multiproc = False
            if request.POST.get("multiproc") != None:
                multiproc = True
            

            from .tasks import run_classification

            
            current_user = request.user
            current_user_email = current_user.email
            current_user_isauth = current_user.is_authenticated
            
            if request.method in ('GET','POST'):
                run_classification.delay(dct, mstid_index, new_list, 
                                            recompute, reupdate_db, music_process, 
                                            music_new_list, music_reupdate_db,
                                            nprocs, multiproc,current_user_email,current_user_isauth)
            # from time import sleep
            # sleep(2)
            return redirect("/manual")
    except:
        pass
        
    return render(request,"classify.html")

def manual_search(request):
    mstid_list  = get_active_list()
    # mstid_list  = mongo_tools.get_active_list()

    mstidDayDict,quietDayDict,noneDayDict,unclassifiedDayDict = loadDayLists(mstid_list=mstid_list)
    # mstidDayDict,quietDayDict = mongo_tools.get_mstid_days()
    
    mstidStr  = linkUp(mstidDayDict)
    quietStr  = linkUp(quietDayDict)
    noneStr   = linkUp(noneDayDict)
    unclassifiedStr   = linkUp(unclassifiedDayDict)
    # import ipdb;ipdb.set_trace()
    webData = {}
    webData['mstidDays']      = mstidStr
    webData['quietDays']      = quietStr
    webData['noneDays']       = noneStr
    webData['unclassifiedDays'] = unclassifiedStr 

    # Count number of events. ######################################################
    webData['mstid_days_total'] = len(mstidDayDict)
    mstid_days_manual_checked = 0
    for event in mstidDayDict:
      if 'category_manu' in event:
          mstid_days_manual_checked = mstid_days_manual_checked + 1
    webData['mstid_days_manual_checked'] = mstid_days_manual_checked

    webData['quiet_days_total'] = len(quietDayDict)
    quiet_days_manual_checked = 0
    for event in quietDayDict:
      if 'category_manu' in event:
          quiet_days_manual_checked = quiet_days_manual_checked + 1
    webData['quiet_days_manual_checked'] = quiet_days_manual_checked

    webData['none_days_total'] = len(noneDayDict)
    none_days_manual_checked = 0
    for event in noneDayDict:
      if 'category_manu' in event:
          none_days_manual_checked = none_days_manual_checked + 1
    webData['none_days_manual_checked'] = none_days_manual_checked

    webData['unclassified_days_total'] = len(unclassifiedDayDict)
    unclassified_days_manual_checked = 0
    for event in unclassifiedDayDict:
      if 'category_manu' in event:
          unclassified_days_manual_checked = unclassified_days_manual_checked + 1
    webData['unclassified_days_manual_checked'] = unclassified_days_manual_checked

    webData['days_total']       = webData['mstid_days_total'] + webData['quiet_days_total'] 
    """+ webData['none_days_total'] + webData['unclassified_days_total']"""
    ################################################################################

    webData['list_dropdown']    = listDropDown()
    webData['mstid_list']       = mstid_list
    webData['homeURL']          = '/'

    timestamp=datetime.datetime.utcnow().strftime('%Y%m%d%H%M%S')
    base_url =  "{0}://{1}".format(request.scheme, request.get_host())
    return render(request,'manual.html',{'webData':webData,'timestamp':timestamp,'base_url':base_url})
    

#User home Page
def home(request):
    return render(request,'home.html')



#Plot range time plot
def plot_rtp(request):
    radar   = 'bpk'
    sDate   = datetime.datetime(2017,1,15,1)
    eDate   = datetime.datetime(2017,1,15,23)
    fit_sfx = "fitacf"
    data_dir = f'/home/fran/code/SuperdarnW3usr/ForGitRepo/'
    fitacf  = load_fitacf(radar,sDate,eDate,data_dir=data_dir)
    dataObj = music.musicArray(fitacf,sTime=sDate,eTime=eDate,fovModel='GS')
    fig = plt.figure(figsize=(20,10))
    ax  = fig.add_subplot(121)
    img = BytesIO()
    musicRTP(dataObj,axis=ax,beam=13)
    plt.savefig(img,format='png')
    plt.close()
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode('utf8')
    img.close()
    return render(request,'test.html',{'plot_url':plot_url})

#Do manual search

# Works

def select_source(request):
    result = 1
    new_source    = request.GET['sources_dropdown']
    source_link   = msc.get_output_path()
    if os.path.islink(source_link):
        os.remove(source_link)
        os.symlink(new_source,source_link)
  
        flash('Data source selected: {}'.format(new_source))
        result=0
    
    return JsonResponse(result=result)

# Works
# @app.route('/list_save_as',methods=['GET'])
def list_save_as(request):
    list_name   = request.GET['listName']
    saveAsName  = request.GET['saveAsName']
    mstid_list  = request.GET['mstid_list']
       
    if list_name == 'saveAs':
        saveAs    = True
        list_name  = saveAsName
    else:
        saveAs    = False

    checkName = list_name in db.list_collection_names()
    if list_name == '' or list_name == 'clean':
        result = 1
    elif checkName == True and saveAs == True:
        flash('Error! A collection with the name "'+list_name+'" already exists!')
        result = 1
    else:
        if checkName == False:
            db['listTracker'].insert_one({'name':list_name})
            for x in db[mstid_list].find():
                db[list_name].insert_one(x)
            set_active_list(list_name)
            result = 0
    result_state_list= {}
    result_state_list['result'] = result
    return JsonResponse(result_state_list)

# Works
# @app.route('/load_list', methods=['GET'])
def load_list(request):
  listName    = request.GET['list_dropdown']
  mstid_list  = request.GET['mstid_list']
  result = 1
  if listName != '' and listName != 'saveAs':
    set_active_list(listName)
    # flash('List "'+listName+'" loaded.')
    messages.add_message(request, messages.SUCCESS, 'List {0} loaded'.format(listName))
    result=0
  result_state_load            = {}
  result_state_load['result'] = 0
  return JsonResponse(result_state_load)

# Works
# @app.route('/list_delete', methods=['GET'])
def list_delete(request):
    listName    = request.GET['list_dropdown']
    mstid_list  = request.GET['mstid_list']
    result = 1
    if listName != '' and listName != 'saveAs' and listName != None and listName != 'clean':
        entry = db['listTracker'].delete_one({"name":listName})
        db[listName].drop()
        # flash('List "'+listName+'" deleted.')
        messages.add_message(request, messages.SUCCESS, 'List {0} deleted'.format(listName))
        result=0
    result_state_delete = {}
    result_state_delete['result'] = result
    return JsonResponse(result_state_delete)
#works
# @app.route('/update_category',methods=['GET'])
def update_category(request):
    '''Update the categories datebase for a particular day.'''

    #Unpack variables from the post.
    mstid_list                = request.GET['mstid_list']
    settings_collection_name  = request.GET['settings_collection_name']

    item                  = {}
    item['radar']         = request.GET['categ_radar']
    item['date']          = datetime.datetime.strptime(request.GET['categ_day'],'%Y%m%d-%H%M')

    category_manu = request.GET['categ_manu']
    if category_manu != 'Null':
        if category_manu != 'None':
            item['category_manu'] = category_manu
        else:
            item['category_manu'] = 'None'

    checked = request.GET['categ_checked']
    if checked == 'true':
        item['checked'] = True
    else:
        item['checked'] = False

#    item['notes']         = request.GET['categ_notes', 0, type=str)

    tmp = db_update_mstid_list(item,mstid_list=mstid_list)
    # tmp = updateDb_mstid_list()
    radar = item['radar']
    mstidDayDict,quietDayDict,noneDayDict,unclassifiedDayDict = loadDayLists(mstid_list=mstid_list)

    mstidStr  = linkUp(mstidDayDict)
    quietStr  = linkUp(quietDayDict)
    noneStr   = linkUp(noneDayDict)
    unclassifiedStr   = linkUp(unclassifiedDayDict)
    result_mstid_Str = {}
    result_mstid_Str['mstidStr'] = mstidStr
    result_mstid_Str['quietStr'] = quietStr
    result_mstid_Str['noneStr']  = noneStr

    return JsonResponse(result_mstid_Str)
#HOMEPAGE ends

# # @app.route('/music_update_category',methods=['GET'])
# def music_update_category(request):
#     '''Update the categories datebase for a particular day.
#        This version is called from the music_edit page.'''

#     #Unpack variables from the post.
#     mstid_list      = request.GET['mstid_list']
#     category_manu   = request.GET['categ_manu']
#     str_id          = request.GET['_id']

#     _id = ObjectId(str_id)
#     event   = db[mstid_list].find_one({'_id':_id})

#     status = db[mstid_list].update_one({'_id':_id},{'$set': {'category_manu':category_manu}})
#     # import ipdb;ipdb.set_trace()

#     if 'err' not in status:
#         result=0
#     else:
#         result=1

#     nav_mode = get_nav_mode()
#     prev_url,next_url  = msc.get_prev_next(mstid_list,_id,mode=nav_mode)
#     urls= {}
#     urls['prev_url'] = prev_url
#     urls['next_url'] = next_url
#     music_update_category_result = {}
#     music_update_category_result['result'] = result
#     music_update_category_result['category_manu'] = category_manu
#     music_update_category_result['urls'] = urls
#     return JsonResponse(music_update_category_result)

# @app.route('/update_nav_mode',methods=['GET'])
def update_nav_mode(request):

    #Unpack variables from the post.
    nav_mode        = request.GET['nav_mode','list']
    mstid_list      = request.GET['mstid_list']
    str_id          = request.GET['_id']
    _id = ObjectId(str_id)

    set_nav_mode(nav_mode)

    urls = msc.get_prev_next(mstid_list,_id,mode=nav_mode)
    update_nav_mode_result = {}
    update_nav_mode_result['prev_url'] = urls[0]
    update_nav_mode_result['next_url'] = urls[1]
    
    return JsonResponse(update_nav_mode_result)


# Generate RTI Plot ############################################################
def plot_rti(request):
  '''Plot an RTI plot for the given day to a PNG and return the PNG's location and
  some information about that day from the database.'''

  #Unpack variables from the get.
  radar       = request.GET['radar']
  gwDay       = request.GET['gwDay']
  param       = request.GET['param']
  mstid_list  = request.GET['mstid_list']
  
  stime   = datetime.datetime.strptime(gwDay,'%Y%m%d-%H%M')
  etime   = stime + datetime.timedelta(hours=2)
  shortDt = datetime.datetime.strptime(gwDay[0:8],'%Y%m%d')
  
  #Build the path of the RTI plot and call the plotting routine.
  d = '/staticfiles/rti'
  outputFile = d+'/'+stime.strftime('%Y%m%d')+'.'+radar+'.'+param+'rti.png'

  try:
    os.makedirs(d)
  except:
    pass
  # import ipdb;ipdb.set_trace()

  if not os.path.exists(outputFile):
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 22}
    matplotlib.rc('font', **font)
    # tick_size   = 16

    # xticks  = []
    # hours   = np.arange(0,26,2) #Put a tick every 2 hours.
    # for hour in hours:
    #     tmp = shortDt + datetime.timedelta(hours=float(hour))
    #     xticks.append(tmp)
    # axvlines  = xticks[:]
    from pyDARNmusic import load_fitacf
    from pyDARNmusic.music import musicArray
    fitacf  = load_fitacf(radar,stime,etime)
    dataObj = musicArray(fitacf,sTime=stime,eTime=etime,fovModel='GS')
    fig = Figure(figsize=(20,10))
    ax  = fig.add_subplot(111)
    musicRTP(dataObj,axis=ax)
    # fig = musicRTP(shortDt,radar,params=['power'],show=False,retfig=True,figure=fig,xtick_size=tick_size,ytick_size=tick_size,xticks=xticks,axvlines=axvlines)
    # canvas = FigureCanvasAgg(fig)
    
    # fig.savefig("webserver"+outputFile)
    fig.savefig(outputFile[1:])

  #Load in infomation about day stored in database.
  dbDict = db[mstid_list].find_one({'radar':radar,'date':stime})
  if dbDict == None: dbDict = {}

  #Pack everything up into a dictionary to return to the JavaScript function.
  output            = {}
  # import ipdb;ipdb.set_trace()
  #Append timestamp to force reload of image.
  output['result']  = outputFile+datetime.datetime.now().strftime('?%Y%m%d%H%M%S')
#   output['result'] = d+"?"+ stime.strftime('%Y%m%d')+'.'+radar+'.'+param+'rti.png'
  output['radar']   = radar
  output['gwDay']   = gwDay

  output['categ_auto_mstid'] = False
  output['categ_auto_quiet'] = False
  output['categ_auto_none']  = False
  if 'category_auto' in dbDict:
    if dbDict['category_auto'] == 'mstid':
      output['categ_auto_mstid'] = True
    elif dbDict['category_auto'] == 'quiet':
      output['categ_auto_quiet'] = True
    elif dbDict['category_auto'] == 'None':
      output['categ_auto_none'] = True

  output['categ_manu_mstid'] = False
  output['categ_manu_quiet'] = False
  output['categ_manu_none']  = False
  if 'category_manu' in dbDict:
    if dbDict['category_manu'] == 'mstid':
      output['categ_manu_mstid'] = True
    elif dbDict['category_manu'] == 'quiet':
      output['categ_manu_quiet'] = True
    elif dbDict['category_manu'] == 'None':
      output['categ_manu_none'] = True

  if 'checked' in dbDict:
    output['categ_checked'] = dbDict['checked']
  else:
    output['categ_checked'] = None

  if 'survey_code' in dbDict:
    output['survey_code'] = dbDict['survey_code']
  else:
    output['survey_code'] = None

  if 'mlt' in dbDict:
    mlt = dbDict['mlt']
    hr  = int(mlt) * 100.
    mn  = (mlt%1) * 60.
    mlt_str = '%04.0f' % (hr+mn)
    output['mlt'] = mlt_str
  else:
    output['mlt'] = None

  if 'gscat' in dbDict:
    output['gscat'] = dbDict['gscat']
  else:
    output['gscat'] = None

  if 'notes' in dbDict:
    output['categ_notes'] = dbDict['notes']
  else:
    output['categ_notes'] = None
  # import ipdb;ipdb.set_trace()
  return JsonResponse(output)


# Everthing Above here works


def music_update_category(request):
    '''Update the categories datebase for a particular day.
       This version is called from the music_edit page.'''

    #Unpack variables from the post.
    mstid_list      = request.GET['mstid_list']
    category_manu   = request.GET['categ_manu']
    str_id          = request.GET['_id']

    _id = ObjectId(str_id)
    event   = db[mstid_list].find_one({'_id':_id})

    status = db[mstid_list].update_one({'_id':_id},{'$set': {'category_manu':category_manu}})
    # import ipdb;ipdb.set_trace()

    if 'err' not in status.raw_result:
        # messages.add_message(request, messages.SUCCESS, f'Updated Category to {0}'.format(category_manu))
        result=0
    else:
        result=1

    nav_mode = get_nav_mode()
    prev_url,next_url  = msc.get_prev_next(mstid_list,_id,mode=nav_mode)
    music_update_category_urls= {}
    music_update_category_urls['result'] = result
    music_update_category_urls['category_manu'] = category_manu
    music_update_category_urls['prev_url'] = prev_url
    music_update_category_urls['next_url'] = next_url
    # import ipdb;ipdb.set_trace()
    return JsonResponse(music_update_category_urls)


def update_nav_mode(request):

    #Unpack variables from the post.
    nav_mode        = request.GET['nav_mode']
    if not nav_mode:
      nav_mode = request.GET['list']
    mstid_list      = request.GET['mstid_list']
    str_id          = request.GET['_id']
    _id = ObjectId(str_id)

    set_nav_mode(nav_mode)

    urls = msc.get_prev_next(mstid_list,_id,mode=nav_mode)
    result_update_mode = {}
    result_update_mode['prev_url'] = urls[0]
    result_update_mode['next_url'] = urls[1]
    
    return JsonResponse(result_update_mode)



def music(request):
    
    mstid_list  = mongo_tools.get_active_list()

    mstidDayDict,quietDayDict,noneDayDict,unclassifiedDayDict = mongo_tools.loadDayLists(mstid_list=mstid_list)

    mstidStr  = msc.linkUp(mstidDayDict)
    quietStr  = msc.linkUp(quietDayDict)
    noneStr   = msc.linkUp(noneDayDict)
    unclassifiedStr   = msc.linkUp(unclassifiedDayDict)

    webData = {}

    webData['mstidDays']      = mstidStr
    webData['quietDays']      = quietStr
    webData['noneDays']       = noneStr
    webData['unclassifiedDays'] = unclassifiedStr 

    # Count number of events. ######################################################
    webData['mstid_days_total'] = len(mstidDayDict)
    mstid_days_manual_checked = 0
    for event in mstidDayDict:
      if 'category_manu' in event:
          mstid_days_manual_checked = mstid_days_manual_checked + 1
    webData['mstid_days_manual_checked'] = mstid_days_manual_checked

    webData['quiet_days_total'] = len(quietDayDict)
    quiet_days_manual_checked = 0
    for event in quietDayDict:
      if 'category_manu' in event:
          quiet_days_manual_checked = quiet_days_manual_checked + 1
    webData['quiet_days_manual_checked'] = quiet_days_manual_checked

    webData['none_days_total'] = len(noneDayDict)
    none_days_manual_checked = 0
    for event in noneDayDict:
      if 'category_manu' in event:
          none_days_manual_checked = none_days_manual_checked + 1
    webData['none_days_manual_checked'] = none_days_manual_checked

    webData['unclassified_days_total'] = len(unclassifiedDayDict)
    unclassified_days_manual_checked = 0
    for event in unclassifiedDayDict:
      if 'category_manu' in event:
          unclassified_days_manual_checked = unclassified_days_manual_checked + 1
    webData['unclassified_days_manual_checked'] = unclassified_days_manual_checked

    webData['days_total']       = webData['mstid_days_total'] + webData['quiet_days_total'] + webData['none_days_total'] + webData['unclassified_days_total']
    ################################################################################

    # Aggregate things to get out important date info. #############################
    totalDayDictList    = mstidDayDict + quietDayDict + noneDayDict + unclassifiedDayDict
    totalMLTList        = [x['mlt'] for x in totalDayDictList]
    mstidMLTList        = [x['mlt'] for x in mstidDayDict]
    totalMLTArr         = np.array(totalMLTList)
    mstidMLTArr         = np.array(mstidMLTList)

    webData['list_dropdown']    = listDropDown()
    webData['mstid_list']       = mstid_list
    webData['homeURL']          = '/'
    
    try:
        webData['min_mlt']  = '{0:.2f}'.format(totalMLTArr.min())
        webData['max_mlt']  = '{0:.2f}'.format(totalMLTArr.max())
    except:
        pass

    webData['nr_mstid_08_10MLT']    = '{0:d}'.format(np.sum(np.logical_and(mstidMLTArr >=  8., mstidMLTArr < 10.)))
    webData['nr_mstid_10_12MLT']    = '{0:d}'.format(np.sum(np.logical_and(mstidMLTArr >= 10., mstidMLTArr < 12.)))
    webData['nr_mstid_12_14MLT']    = '{0:d}'.format(np.sum(np.logical_and(mstidMLTArr >= 12., mstidMLTArr < 14.)))
    webData['nr_mstid_14_16MLT']    = '{0:d}'.format(np.sum(np.logical_and(mstidMLTArr >= 14., mstidMLTArr < 16.)))
    webData['nr_mstid_16_18MLT']    = '{0:d}'.format(np.sum(np.logical_and(mstidMLTArr >= 16., mstidMLTArr < 18.)))

    messages    = []
    if len(messages) > 0:
        webData['messages_on'] = True
    else:
        webData['messages_on'] = False
    webData['messages']         = messages

    webData['source_selector']      = msc.sourcesDropDown()
    
    timestamp=datetime.datetime.utcnow().strftime('%Y%m%d%H%M%S')
    # tmp = url_for('music_edit',name='name')
    base_url =  "{0}://{1}".format(request.scheme, request.get_host())
    # import ipdb;ipdb.set_trace()
    return render(request,'music.html',{'webData':webData,'timestamp':timestamp,'base_url':base_url})


def music_edit(request):
    '''Plot an RTI plot for the given day to a PNG and return the PNG's location and
    some information about that day from the database.'''
    mstid_list  = get_active_list()

    #Define some dictionaries to be sent to the web template,
    params  = {}
    webData = {}
    musicParams = {}

    #Unpack variables from the get.
    radar       = request.GET['radar']
    sDate       = request.GET['sDate']
    fDate       = request.GET['fDate']
    _idStr      = request.GET['id']

    _id         = ObjectId(_idStr)
    rec         = db[mstid_list].find_one({'_id': _id})
    #  param       = request.GET['param', 0, type=str)
    #  mstid_list  = request.GET['mstid_list', None, type=str)


    sDatetime = datetime.datetime.strptime(sDate,'%Y%m%d.%H%M') 
    fDatetime = datetime.datetime.strptime(fDate,'%Y%m%d.%H%M')

    #Build up list of everything in database record to make it easy to print
    #everything out in Jinja2.
    record_list = []
    keys = list(rec.keys())
    keys.sort()
    for key in keys:
        record_list.append({'key':key,'value':rec[key]})

    musicParams_list = False    #Don't show runfile list in web if there is no run file.
    try:
        runFile     = msc.load_runfile(radar,sDatetime,fDatetime)
        musicParams = runFile.runParams

        #List for Jinja2 dump.
        musicParams_list = []
        keys = list(musicParams.keys())
        keys.sort()
        for key in keys:
            musicParams_list.append({'key':key,'value':musicParams[key]})

    except:
        musicParams['radar']                = radar
        musicParams['sDatetime']            = sDatetime
        musicParams['fDatetime']            = fDatetime

        musicParams['beamLimits']           = None
        musicParams['gateLimits']           = None
        musicParams['interpRes']            = 120
        musicParams['filter_numtaps']       = 101
        musicParams['filter_cutoff_low']    = 0.0003
        musicParams['filter_cutoff_high']   = 0.0012
        musicParams['kx_max']               = 0.05
        musicParams['ky_max']               = 0.05
        musicParams['autodetect_threshold'] = 0.35
        musicParams['neighborhood']         = (10, 10)
    
    if 'beamLimits' in musicParams:
        if np.size(musicParams['beamLimits']) == 2:
            musicParams['beamLimits_0'] = musicParams['beamLimits'][0]
            musicParams['beamLimits_1'] = musicParams['beamLimits'][1]
        else:
            musicParams['beamLimits_0'] = None
            musicParams['beamLimits_1'] = None
    else:
        musicParams['beamLimits_0'] = None
        musicParams['beamLimits_1'] = None

    if 'gateLimits' in musicParams:
        if np.size(musicParams['gateLimits']) == 2:
            musicParams['gateLimits_0'] = musicParams['gateLimits'][0]
            musicParams['gateLimits_1'] = musicParams['gateLimits'][1]
        else:
            musicParams['gateLimits_0'] = None
            musicParams['gateLimits_1'] = None
    else:
        musicParams['gateLimits_0'] = None
        musicParams['gateLimits_1'] = None

    if 'neighborhood' in musicParams:
        if np.size(musicParams['neighborhood']) == 2:
            musicParams['neighborhood_0'] = musicParams['neighborhood'][0]
            musicParams['neighborhood_1'] = musicParams['neighborhood'][1]
        else:
            musicParams['neighborhood_0'] = 10
            musicParams['neighborhood_1'] = 10
    else:
        musicParams['neighborhood_0'] = 10
        musicParams['neighborhood_1'] = 10

    sTime   = musicParams.get('sDatetime',musicParams.get('sTime'))
    eTime   = musicParams.get('fDatetime',musicParams.get('eTime'))

    musicPath   = msc.get_output_path(musicParams['radar'],sTime,eTime)
    pickleName  = msc.get_pickle_name(musicParams['radar'],sTime,eTime)
    picklePath  = os.path.join(musicPath,pickleName)

    dataObj = None
    try:
        with open(picklePath,'rb') as fl:
            dataObj     = pickle.load(fl)
            webData['musicObjStatusClass']  = 'statusNormal'
            webData['musicObjStatus']       = 'Using musicObj file < '+ picklePath +' >.'
    except:
            webData['musicObjStatusClass']  = 'warning'
            webData['musicObjStatus']       = 'MusicObj does not exist: %s' % picklePath

    no_data = False
    if hasattr(dataObj,'messages'):
        if 'No data for this time period.' in dataObj.messages:
            no_data = True
            webData['good_period_warn']       = 'No data for time period. (%s)' % picklePath

    if webData['musicObjStatusClass'] == 'statusNormal' and not no_data:
        dataObj     = checkDataQuality(dataObj,dataSet='originalFit',sTime=sDatetime,eTime=fDatetime)
        dataSets    = dataObj.get_data_sets()
        lst = []
        for dataSet in dataSets:
            currentData = getattr(dataObj,dataSet)
            ds = {}
            ds['name']    = dataSet
            
            histList = []
            keys = list(currentData.history.keys())
            keys.sort()
            for key in keys:
                histList.append({'name':key,'value':currentData.history[key]})
            ds['history'] = histList

            metaList = []
            keys = list(currentData.metadata.keys())
            keys.sort()
            for key in keys:
                metaList.append({'name':key,'value':currentData.metadata[key]})
            ds['metadata'] = metaList

            lst.append(ds)
        webData['dataSets'] = lst

        #Send information about detected signals to the web.
        currentData = getattr(dataObj,dataSets[-1])
        if hasattr(currentData,'sigDetect'):
            sigs        = currentData.sigDetect
            webData['sigList'] = (sigs.string())

        #Stringify information about signals aready in the database...
        if 'signals' in rec:
            webData['sigsInDb'] = stringify_signal_list(rec['signals'],sort_key='serialNr')

        #Tell web if marked as a bad period.
        try:
            if not dataObj.DS000_originalFit.metadata['good_period']:
                webData['good_period_warn'] = 'WARNING: Data marked as bad period!!'
        except:
            pass


    webData['rtiplot_sDatetime'], webData['rtiplot_fDatetime']  = msc.get_default_rti_times(musicParams,dataObj)
    webData['rtiplot_yrange0'], webData['rtiplot_yrange1']      = msc.get_default_gate_range(musicParams,dataObj)
    rti_beam_list = msc.get_default_beams(musicParams,dataObj)

    rti_beam_list_str = [str(rtibm) for rtibm in rti_beam_list]
    webData['rtiplot_beams'] = ','.join(rti_beam_list_str)

    if webData['rtiplot_yrange0'] is None: webData['rtiplot_yrange0'] = 'None'
    if webData['rtiplot_yrange0'] is None: webData['rtiplot_yrange1'] = 'None'

    #See if RTI Plot exists... if so, show it!
    rtiPath     = os.path.join(musicPath,'000_originalFit_RTI.png')
    rtiStaticUrl = "/"+rtiPath
    try:
        with open(rtiPath):
            webData['rtiPath']  = rtiPath[10:] #save rtp url without webserver/
            webData["rtiStaticURL"] = rtiStaticUrl
    except:
        pass

    #If kArr.png exists, show it on top!
    karrPath    = glob.glob(os.path.join(musicPath,'*karr.png'))
    if len(karrPath) > 0:
        webData['karrPath'] = karrPath[0][10:]
    else:
        pass

    #Show all other plots...
    plots = glob.glob(os.path.join(musicPath,'*.png'))
    if plots == []:
        plotDictList = False 
    else:
        plots.sort()
        plotDictList = []
        for plot in plots:
            plotDict = {}
            plotDict['path'] = "/"+plot #removes the webserver/ in path
            plotDict['basename'] = os.path.basename(plot)
            plotDictList.append(plotDict)
        
    webData['plots']        = plotDictList
    webData['radar']        = radar
    webData['sDate']        = sDate
    webData['fDate']        = fDate
    webData['mstid_list']   = mstid_list
    webData['homeURL']      = '/music'

    #Send the categ_manu information in webData since that info may not be included in the record info.
    webData['categ_manu_mstid'] = ''
    webData['categ_manu_quiet'] = ''
    webData['categ_manu_none']  = ''
    webData['categ_manu']       = 'Null'

    if 'category_manu' in rec:
        webData['categ_manu'] = rec['category_manu']
        if rec['category_manu'] == 'mstid': webData['categ_manu_mstid']  = 'checked'
        if rec['category_manu'] == 'quiet': webData['categ_manu_quiet']  = 'checked'
        if rec['category_manu'] == 'None':  webData['categ_manu_none']   = 'checked'

    #Computer prev/next urls
    nav_mode = get_nav_mode()
    if nav_mode == 'list':
        webData['nav_mode_list'] = 'checked'
    else:
        webData['nav_mode_list'] = ''

    if nav_mode == 'category':
        webData['nav_mode_categ'] = 'checked'
    else:
        webData['nav_mode_categ'] = ''

    urls = msc.get_prev_next(mstid_list,_id,mode=nav_mode)
    webData['prev_url'] = urls[0]
    webData['next_url'] = urls[1]

    webData['event_dir_url'] = 'http://sd-work1.ece.vt.edu/data/mstid/statistics/webserver/'+musicPath
    webData['source_selector']      = msc.sourcesDropDown()
    enabled_sources = msc.get_enabled_sources()

    timestamp=datetime.datetime.utcnow().strftime('%Y%m%d%H%M%S')
    recid=rec["_id"]
    del rec["_id"]
    rec["id"]=recid
    # import ipdb;ipdb.set_trace()
    return render(request,'music_edit.html'
            ,{'webData':webData
            ,'params':params
            ,'timestamp':timestamp
            ,'record_list':record_list
            ,'record':rec
            ,'musicParams':musicParams
            ,'musicParams_list':musicParams_list})


def create_music_obj(request):
    radar                   = request.GET['radar']
    sTime                   = request.GET['sTime']
    eTime                   = request.GET['eTime']
    beamLimits_0            = request.GET['beamLimits_0']
    beamLimits_1            = request.GET['beamLimits_1']
    gateLimits_0            = request.GET['gateLimits_0']
    gateLimits_1            = request.GET['gateLimits_1']
    interpolationResolution = request.GET['interpolationResolution']
    filterNumtaps           = request.GET['filterNumtaps']
    firFilterLimits_0       = request.GET['firFilterLimits_0']
    firFilterLimits_1       = request.GET['firFilterLimits_1']
    window_data             = request.GET['window_data']
    kx_max                  = request.GET['kx_max']
    ky_max                  = request.GET['ky_max']
    autodetect_threshold_str = request.GET['autodetect_threshold']
    neighborhood_0          = request.GET['neighborhood_0']
    neighborhood_1          = request.GET['neighborhood_1']

    #Convert string type before sending to music object creation.
    sDatetime = datetime.datetime.strptime(sTime,'%Y-%m-%d %H:%M:%S')
    fDatetime = datetime.datetime.strptime(eTime,'%Y-%m-%d %H:%M:%S')

    current_user = request.user
    current_user_email = current_user.email
    current_user_isauth = current_user.is_authenticated
    from .tasks import create_music_object_task
    create_music_object_task(radar, sDatetime, fDatetime,
                             beamLimits_0, beamLimits_1, gateLimits_0,
                             gateLimits_1, interpolationResolution,
                             filterNumtaps, firFilterLimits_0, 
                             firFilterLimits_1, window_data,
                             kx_max,ky_max, autodetect_threshold_str,
                             neighborhood_0, neighborhood_1,
                             current_user_email, current_user_isauth)

    # from django.core.management import call_command
    # import ipdb;ipdb.set_trace()
    # call_command("collectstatic", interactive=False)

    # music_obj_result={}
    # music_obj_result['result'] = 0
    result = {}
    result["result"]=0
    return JsonResponse(result)


def run_music(request):
    runfile_path    = request.GET['runfile_path']
    # import ipdb;ipdb.set_trace()
    # msc.run_music(runfile_path)
    # msc.music_plot_all(runfile_path)
    current_user = request.user
    current_user_email = current_user.email
    current_user_isauth = current_user.is_authenticated
    
    from .tasks import run_music_task
    run_music_task(runfile_path, current_user_email, current_user_isauth)
    
    run_music_result = {}
    run_music_result['result'] = 0
    
    return JsonResponse(run_music_result)
   
    # return render(request, '404.html', status=404)


def music_plot_all(request):
    runfile_path    = request.GET['runfile_path']
    msc.music_plot_all(runfile_path)

    music_plot_all_result = {}
    music_plot_all_result['result'] = 0
    return JsonResponse(music_plot_all_result)


def music_plot_fan(request):
    runfile_path    = request.GET['runfile_path']
    mstid_list      = request.GET['mstid_list']
    str_id          = request.GET['_id']
    time            = request.GET['time']
    sDatetime       = datetime.datetime.strptime(time,'%Y-%m-%d %H:%M:%S')
    fanScale_0            = request.GET['fanScale_0']
    fanScale_1            = request.GET['fanScale_1']

    try:
        fanScale = (float(fanScale_0), float(fanScale_1))
    except:
        fanScale = None

    fileName        = '001_beamInterp_fan.png'
    msc.music_plot_fan(runfile_path,time=sDatetime,fileName=fileName,scale=fanScale)

    music_plot_fan_result = {}
    music_plot_fan_result['result'] = 0
    return JsonResponse(music_plot_fan_result)
#WORK!!!!

def music_plot_rti(request):
    runfile_path    = request.GET['runfile_path']
    mstid_list      = request.GET['mstid_list']
    str_id          = request.GET['_id']
    beam            = request.GET['beam']
    sTime           = request.GET['sTime']
    eTime           = request.GET['eTime']
    sDatetime       = datetime.datetime.strptime(sTime,'%Y-%m-%d %H:%M:%S')
    eDatetime       = datetime.datetime.strptime(eTime,'%Y-%m-%d %H:%M:%S')
    rtiScale_0      = request.GET['rtiScale_0']
    rtiScale_1      = request.GET['rtiScale_1']
    rtiYrange_0     = request.GET['rtiYrange_0']
    rtiYrange_1     = request.GET['rtiYrange_1']

    beam            = beam.split(',')
    beam            = [int(xx) for xx in beam]

    try:
        rtiScale = (float(rtiScale_0), float(rtiScale_1))
    except:
        rtiScale = None

    try:
        rtiYrange = (float(rtiYrange_0), float(rtiYrange_1))
    except:
        rtiYrange = None

    runFile         = msc.load_runfile_path(runfile_path)
    musicParams     = runFile.runParams
    musicObj_path   = musicParams['musicObj_path']
    dataObj         = pickle.load(open(musicObj_path,'rb'))

    xlim = (sDatetime,eDatetime)

    #Mark where sampling window starts and stops.
    dataObj.DS000_originalFit.metadata['timeLimits'] = [musicParams['sDatetime'],musicParams['fDatetime']]

    rtiPath     = os.path.join(musicParams['path'],'000_originalFit_RTI.png')
    msc.plot_music_rti(dataObj,fileName=rtiPath,dataSet="originalFit",beam=beam,xlim=xlim,ylim=rtiYrange,scale=rtiScale)

    music_plot_rti_result = {}
    music_plot_rti_result['result'] = 0
    return JsonResponse(music_plot_rti_result)


def add_music_params_db(request):
    #Get data from the webpage.
    runfile_path    = request.GET['runfile_path']
    mstid_list      = request.GET['mstid_list']
    str_id          = request.GET['_id']
    signals         = request.GET['signals']

    #Parse list of signals.
    signal_order_list = [int(x) for x in signals.split(',')]
    signal_order_list.sort()

    #Load the runfile and the associated musicObj.
    runfile         = msc.load_runfile_path(runfile_path)
    picklePath      = runfile.runParams['musicObj_path']

    dataObj     = pickle.load(open(picklePath,'rb'))
    dataSets    = dataObj.get_data_sets()
    currentData = getattr(dataObj,dataSets[-1])

    #Pull up database record to see if there are already items stored.
    _id = ObjectId(str_id)
    event   = db[mstid_list].find_one({'_id':_id})

    if 'signals' in event:
        sigList     = event['signals']
        try:
            serialNr    = max([x['serialNr'] for x in sigList]) + 1
        except:
            serialNr    = 0
    else:
        sigList = []
        serialNr    = 0

#    serialNr = 0
#    sigList = []
    if hasattr(currentData,'sigDetect'):
        sigs    = currentData.sigDetect
        for order in signal_order_list:
            for sig in sigs.info:
                if sig['order'] == order:
                    sigInfo = {}
                    sigInfo['order']    = int(sig['order'])
                    sigInfo['kx']       = float(sig['kx'])
                    sigInfo['ky']       = float(sig['ky'])
                    sigInfo['k']        = float(sig['k'])
                    sigInfo['lambda']   = float(sig['lambda'])
                    sigInfo['azm']      = float(sig['azm'])
                    sigInfo['freq']     = float(sig['freq'])
                    sigInfo['period']   = float(sig['period'])
                    sigInfo['vel']      = float(sig['vel'])
                    sigInfo['max']      = float(sig['max'])
                    sigInfo['area']     = float(sig['area'])
                    sigInfo['serialNr'] = serialNr
                    sigList.append(sigInfo)
                    serialNr = serialNr + 1

    status = db[mstid_list].update_one({'_id':_id},{'$set': {'signals':sigList}})

    add_music_params_db_result = {}
    add_music_params_db_result['result'] = 0
    return JsonResponse(add_music_params_db_result)


def del_music_params_db(request):
    #Get data from the webpage.
    runfile_path    = request.GET['runfile_path']
    mstid_list      = request.GET['mstid_list']
    str_id          = request.GET['_id']
    signals         = request.GET['signals']

    #Parse list of signals.
    signal_serialNr_list = [int(x) for x in signals.split(',')]

    #Pull up database record to see if there are already items stored.
    _id     = ObjectId(str_id)
    event   = db[mstid_list].find_one({'_id':_id})
    sigList = event['signals']

    for sig in list(sigList):
        if sig['serialNr'] in signal_serialNr_list:
            sigList.remove(sig)

    status  = db[mstid_list].update_one({'_id':_id},{'$set': {'signals':sigList}})

    del_music_params_db_result = {}
    del_music_params_db_result['result'] = 0
    return JsonResponse(del_music_params_db_result)



def update_music_analysis_status(request):
    #Get data from the webpage.
    runfile_path    = request.GET['runfile_path']
    mstid_list      = request.GET['mstid_list']
    str_id          = request.GET['_id']
    analysis_status = request.GET['analysis_status']

    #Pull up database record to see if there are already items stored.
    _id     = ObjectId(str_id)
    event   = db[mstid_list].find_one({'_id':_id})

    status  = db[mstid_list].update_one({'_id':_id},{'$set': {'music_analysis_status':bool(analysis_status)}})
    update_music_analysis_status_result = {}
    update_music_analysis_status_result['result'] = 0
    return JsonResponse(update_music_analysis_status_result)


def add_to_detected(request):
    #Get data from the webpage.
    runfile_path    = request.GET['runfile_path']
    mstid_list      = request.GET['mstid_list']
    str_id          = request.GET['_id']
    new_kx          = request.GET['new_kx']
    new_ky          = request.GET['new_ky']

    if new_kx == None: 
        return jsonify(result=0)
    if new_ky == None: 
        return jsonify(result=0)

    #Load the runfile and the associated musicObj.
    runfile         = msc.load_runfile_path(runfile_path)
    picklePath      = runfile.runParams['musicObj_path']
    musicPath       = runfile.runParams['path']

    dataObj     = pickle.load(open(picklePath,'rb'))
    add_signal(new_kx,new_ky,dataObj,dataSet='active')
    pickle.dump(dataObj,open(picklePath,'wb'))

    karrPath    = glob.glob(os.path.join(musicPath,'*karr.png'))[0]
    msc.music_plot_karr(runfile_path,karrPath)

    add_to_detected_result = {}
    add_to_detected_result['result'] = 0
    return JsonResponse(add_to_detected_result)


def del_from_detected(request):
    #Get data from the webpage.
    runfile_path    = request.GET['runfile_path']
    mstid_list      = request.GET['mstid_list']
    str_id          = request.GET['_id']
    signals         = request.GET['signals']

    #Parse list of signals.
    signal_order_list = [int(x) for x in signals.split(',')]
    signal_order_list.sort()

    #Load the runfile and the associated musicObj.
    runfile         = msc.load_runfile_path(runfile_path)
    picklePath      = runfile.runParams['musicObj_path']
    musicPath       = runfile.runParams['path']

    dataObj     = pickle.load(open(picklePath,'rb'))
    del_signal(signal_order_list,dataObj,dataSet='active')
    pickle.dump(dataObj,open(picklePath,'wb'))

    karrPath    = glob.glob(os.path.join(musicPath,'*karr.png'))[0]
    msc.music_plot_karr(runfile_path,karrPath)

    del_from_detected_result = {}
    del_from_detected_result['result'] = 0
    return JsonResponse(del_from_detected_result)

def page_not_found_view(request, exception):
    return render(request, '404.html', status=404)