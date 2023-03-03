from django.shortcuts import render, redirect
from django.contrib.auth import login,logout,authenticate
from .forms import RegisterForm

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