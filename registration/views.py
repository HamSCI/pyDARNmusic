from django.shortcuts import render, redirect
from django.contrib.auth import login,logout,authenticate
from django.contrib.auth.forms import UserChangeForm
from .forms import RegisterForm, ProfileForm

# Create your views here.
def register(request):
    if request.method == "POST":
        form = RegisterForm(request.POST)
        if form.is_valid():
            user = form.save()
            login(request,user)
            return redirect("/home")
    else:
        form = RegisterForm()
    return render(request,'register.html',{"form":form})

# def login(request):
def usereditview(request):
    if request.method == "POST":
        form2 = ProfileForm(request.POST, instance=request.user)
        if form2.is_valid():
            form2.save()
            return redirect("/home")
    else:
        form2 = ProfileForm()
    return render(request,'edit_profile.html',{"form2":form2})