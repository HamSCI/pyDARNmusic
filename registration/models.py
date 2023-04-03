from django.db import models
from django.contrib.auth.models import AbstractUser
# from django import forms
# Create your models here.
class UserProfile(AbstractUser):
    USER_TYPE_ATMOSPHERIC_RESEARCH_STUDENT = "Atmospheric Research Student"
    USER_TYPE_ATMOSPHERIC_RESEARCH_SCIENTIST = "Atmospheric Research Scientist"
    USER_TYPE_US_CITIZEN = "US Citizen"
    USER_TYPES = [
        (USER_TYPE_ATMOSPHERIC_RESEARCH_STUDENT, "ATMOSPHERIC RESEARCH STUDENT"),
        (USER_TYPE_ATMOSPHERIC_RESEARCH_SCIENTIST, "ATMOSPHERIC RESEARCH SCIENTIST"),
        (USER_TYPE_US_CITIZEN, "US CITIZEN"),
    ]
    # username = forms.CharField(max_length=250,required=True)
    # email = forms.EmailField(required=True, unique=True)
    # password1 = 
    user_type = models.CharField(max_length=250, choices=USER_TYPES)

class Profile(models.Model):
    user = models.OneToOneField(UserProfile, null=True, on_delete=models.CASCADE)

    def __str__(self):
        return str(self.user)