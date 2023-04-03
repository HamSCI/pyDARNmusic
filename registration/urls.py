from django.urls import path
from django.contrib.auth import views as auth_views
from . import views

urlpatterns = [
    path('register/',views.register, name='register'),
    path('edit_profile/',views.usereditview, name='edit_profile'),
    path('password/',auth_views.PasswordChangeView.as_view(template_name='registration/change-password.html'), name='change_password'),
]
