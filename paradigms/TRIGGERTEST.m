for i = 1:255
sendtrigger(i,port,SITE,stayup)
pause(0.05)
end

for i = 0:7
    val = 2^i
sendtrigger(val,port,SITE,stayup)
pause(0.5)
end

patterns = [1, 3, 7, 15, 31, 63, 127, 255]
for i = length(patterns)
sendtrigger(patterns(i),port,SITE,stayup)
pause(0.5)
end