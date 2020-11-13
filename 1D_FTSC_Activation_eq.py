from si_cfd import SiData
from si_cfd import CalDiff
import matplotlib.pyplot as plt

def plot_data(x, qn, time):
	plt.plot(x, qn, marker='o', lw=2, label=f'time={time}')

plt.figure(figsize=(7, 7), dpi=100)
plt.rcParams["font.size"]=10

data=[]
x=[]
start_x=0
end_x=2.0
dx = 0.02
dt = 0.01

i=0
while start_x + dx*i <= end_x :
	x.append(start_x+dx*i)
	if x[i] <= 1.0 :
		data.append(1)
	else :
		data.append(0)
	i=i+1

time = 0
end_time=0.3
q=SiData()
cal=CalDiff()
q.set_data(dt, dx, data)
plot_data(x, q.data, time)
n=0
plot_dn=2*(0.05/dt)
while time <= end_time:
	n = n + 1
	time = round(n*dt,1)
	#qn = (cal.Dt("forward_diff", q) + cal.Dx("centor_diff",q,1)).solve()
	#qn = (cal.Dt("forward_diff", q) + cal.Dx("advance_diff", q, -1)).solve()
	qn = ( cal.Dt("Lax_diff", q) + cal.Dx("Lax_diff", q, -1)).solve()
	#qn = ( cal.Dt("Lax_Wendroff_diff",q) + cal.Dx("Lax_Wendroff_diff", q, -1) ).solve()
	q.update(qn,time)
	q.boundary()
	if n % plot_dn == 0 :
		plot_data(x, q.data, time)
	#time = round(n*dt,1)
	
#図の後処理
plt.xlim([0, 2])
#plt.ylim([0,1.4])
plt.xlabel('x')
plt.ylabel('q')
plt.legend(loc="bottom left")
plt.show()
