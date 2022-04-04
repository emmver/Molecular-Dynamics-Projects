t=1e-5
req_t=5e-3
for i in range(1000): 
    if t<req_t:
        t+=t*0.2;print(i,t)
    else:
        print('end',i)
        print('t-step',t)
        t=req_t
        print('final t',t)
        break
