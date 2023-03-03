from django.contrib import admin
from django.contrib.auth.admin import UserAdmin
from .models import UserProfile
# Register your models here.
@admin.register(UserProfile)
class MyUserAdmin(UserAdmin):
        model = UserProfile
        list_display = ('username','user_type',
                        'email')
        list_filter = ('username',
                        'email','user_type')
        search_fields = ('username','email' )
        ordering = ('username','email' )