import matplotlib.pyplot as plt

x=[[1,2,3,4],[1,2,3,4],[1,2,3,4],[1,2,3,4]]
y=[[1,2,3,4],[2,3,4,5],[3,4,5,6],[7,8,9,10]]
y2=[[11,12,13,24],[42,33,34,65],[23,54,65,86],[77,90,39,54]]
colours=['r','g','b','k']

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
for i in range(len(x)):
    ax1.plot(x[i],y2[i],colours[i])
    ax2.plot(x[i],y[i],colours[i])

fig1.savefig('fig1.png')
fig2.savefig('fig2.png')
