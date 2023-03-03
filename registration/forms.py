from django.contrib.auth import login, authenticate
from django.contrib.auth.forms import UserCreationForm
# from django.contrib.auth.models import User
from .models import UserProfile
from django import forms


class RegisterForm(UserCreationForm):
    # USER_TYPE_ATMOSPHERIC_RESEARCH_STUDENT = "Atmospheric Research Student"
    # USER_TYPE_ATMOSPHERIC_RESEARCH_SCIENTIST = "Atmospheric Research Scientist"
    # USER_TYPE_US_CITIZEN = "US Citizen"
    # USER_TYPES = [
    #     (USER_TYPE_ATMOSPHERIC_RESEARCH_STUDENT, "ATMOSPHERIC RESEARCH STUDENT"),
    #     (USER_TYPE_ATMOSPHERIC_RESEARCH_SCIENTIST, "ATMOSPHERIC RESEARCH SCIENTIST"),
    #     (USER_TYPE_US_CITIZEN, "US CITIZEN"),
    # ]
    
    # email = forms.EmailField(required=True)
    # user_type = forms.ChoiceField(choices=USER_TYPES,widget=forms.RadioSelect(),required=True)
    
    class Meta:
        model = UserProfile
        fields = ["username","email","user_type","password1","password2"]
        ordering = ['username']