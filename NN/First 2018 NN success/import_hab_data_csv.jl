cd("/home/zachary/Documents/data/Salem HABs/csvs/2.12.19")
using CSV

dates=CSV.read("dates.csv")
tbv1=CSV.read("tbv1.csv")
tbv2=CSV.read("tbv2.csv")
tbv3=CSV.read("tbv3.csv")
tbv4=CSV.read("tbv4.csv")
cyano=CSV.read("tbv5.csv")

temp=CSV.read("temp.csv")
rain=CSV.read("rain.csv")

tox1=CSV.read("tox1.csv")
tox2=CSV.read("tox2.csv")
tox3=CSV.read("tox3.csv")
tox4=CSV.read("tox4.csv")

col=CSV.read("color.csv")
sattime=CSV.read("sattime.csv")
