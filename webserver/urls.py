from django.urls import path
from . import views

urlpatterns = [
    path('',views.apphome),
    path('apphome/',views.apphome, name='apphome'),
    path('home/',views.home,name='userhome'),
    path('plot_rtp/',views.plot_rtp),
    path('manual/',views.manual_search, name='manual'),
    path('classify/',views.classify_mstids, name='classify'),
    path('list_save_as/',views.list_save_as, name='list_save_as'),
    path('select_source/',views.select_source, name='select_source'),
    path('load_list/',views.load_list, name='load_list'),
    path('list_delete/',views.list_delete, name='list_delete'),
    path('update_category/',views.update_category, name='update_category'),
    path('music_update_category/',views.music_update_category, name='music_update_category'),
    path('update_nav_mode/',views.update_nav_mode, name='update_nav_mode'),
    path('rti/',views.plot_rti, name='plot_rti'),
    # MUSIC PAGE
    path('music_update_category/',views.music_update_category, name='music_update_category'),
    path('update_nav_mode/',views.update_nav_mode, name='update_nav_mode'),
    path('music/',views.music, name='music'),
    path('music_edit/',views.music_edit, name='music_edit'),
    path('create_music_obj/',views.create_music_obj, name='create_music_obj'),
    path('run_music/',views.run_music, name='run_music'),
    path('music_plot_all/',views.music_plot_all, name='music_plot_all'),
    path('music_plot_fan/',views.music_plot_fan, name='music_plot_fan'),
    path('music_plot_rti/',views.music_plot_rti, name='music_plot_rti'),
    path('add_music_params_db/',views.add_music_params_db, name='add_music_params_db'),
    path('del_music_params_db/',views.del_music_params_db, name='del_music_params_db'),
    path('update_music_analysis_status/',views.update_music_analysis_status, name='update_music_analysis_status'),
    path('add_to_detected/',views.add_to_detected, name='add_to_detected'),
    path('del_from_detected/',views.del_from_detected, name='del_from_detected'),
    # path('music_edit/static/**')
]
