import numpy as np

class SiSolData():
	An=np.array([]) #求める変数の係数
	Ab=np.zeros([]) #その他の係数
	def __init__(self, dataN):
		self.An=np.zeros((dataN,dataN))
		self.Ab=np.zeros(dataN)
			
	def __add__(self, other):
		self.An = self.An + other.An
		self.Ab = self.Ab + other.Ab
		return self
		
	def solve(self):
		invAn = np.linalg.inv(self.An)
		Ab_T  = np.array([self.Ab]).T
		result = np.dot(invAn, -Ab_T)
		#print ("invAn : \n" ,invAn)
		#print ("Ab : \n", self.Ab)
		#print ("Ab_T : \n" ,-Ab_T)
		#print ("result : \n", np.dot(invAn,-Ab_T))
		return result.T[0]

class CalDiff():
	
	###Dt
	def forward_diff_Dt(self, si_data):
		n=si_data.N
		dt=si_data.dt
		soldata = SiSolData(n)
		soldata.An[0][0]=1
		soldata.Ab[0]=1
		soldata.An[n-1][n-1]=1
		soldata.Ab[n-1]=1
		for i in range(1,n-1):
			soldata.An[i][i]=1.0/dt
			soldata.Ab[i]=(-(1.0/dt)*si_data.data[i])
		return soldata
	
	def Lax_diff_Dt(self, si_data):
		n=si_data.N
		dt=si_data.dt
		soldata = SiSolData(n)
		soldata.An[0][0]=1
		soldata.Ab[0]=1
		soldata.An[n-1][n-1]=1
		soldata.Ab[n-1]=1
		for i in range(1,n-1):
			soldata.An[i][i]=1.0/dt
			soldata.Ab[i]=(-(1.0/dt)*(si_data.data[i+1]+si_data.data[i-1])/2.0)
		return soldata
		
	def Dt(self, sol_name, si_data):
		if sol_name == "forward_diff" :
			return self.forward_diff_Dt(si_data)
		if sol_name == "Lax_diff" :
			return self.Lax_diff_Dt(si_data)
		if sol_name == "Lax_Wendroff_diff":
			return self.forward_diff_Dt(si_data)
		else :
			print("The func type name '{0}' of Dt is not defined.".format(sol_name))
	
	###Dx
	def centor_diff_Dx(self, si_data, coef):
		n=si_data.N
		soldata = SiSolData(n)
		soldata.An[0][0]=0
		soldata.Ab[0]=si_data.data[0]
		soldata.An[n-1][n-1]=0
		soldata.Ab[n-1]=si_data.data[n-1]
		for i in range(1,n-1):
			soldata.Ab[i]=coef*(si_data.data[i+1]-si_data.data[i-1])/(2.0*si_data.dx)
		return soldata
		
	def advance_diff_Dx(self, si_data, coef):
		n=si_data.N
		coef_ab=abs(coef)
		soldata = SiSolData(n)
		soldata.An[0][0]=0
		soldata.Ab[0]=si_data.data[0]
		soldata.An[n-1][n-1]=0
		soldata.Ab[n-1]=si_data.data[n-1]
		for i in range(1,n-1):
			#soldata.Ab[i]=coef*(si_data.data[i]-si_data.data[i-1])/(si_data.dx)
			soldata.Ab[i]=(coef + coef_ab)*0.5*(si_data.data[i]-si_data.data[i-1])/(si_data.dx) + (coef - coef_ab)*0.5*(si_data.data[i+1]-si_data.data[i])/(si_data.dx)
		return soldata
	
	def Lax_Wendroff_diff_Dx(self, si_data, coef):
		n=si_data.N
		dt=si_data.dt
		soldata = SiSolData(n)
		soldata.An[0][0]=0
		soldata.Ab[0]=si_data.data[0]
		soldata.An[n-1][n-1]=0
		soldata.Ab[n-1]=si_data.data[n-1]
		for i in range(1,n-1):
			soldata.Ab[i]=coef*(si_data.data[i+1]-si_data.data[i-1])/(2.0*si_data.dx) - (coef**2)*dt*(si_data.data[i+1]-2.0*si_data.data[i] + si_data.data[i-1])/(2.0*si_data.dx**2)
		return soldata

	#cal.Dx(sol_name, si_data, coef):
	def Dx(self, *args):
		sol_name = args[0]
		si_data = args[1]
		try:
			coef = args[2]
		except:
			coef = 1
			
		if sol_name == "centor_diff" :
			return self.centor_diff_Dx(si_data, coef)
		elif sol_name == "advance_diff" :
			return self.advance_diff_Dx(si_data,coef)
		elif sol_name == "Lax_diff" :
			return self.centor_diff_Dx(si_data, coef)
		elif sol_name == "Lax_Wendroff_diff":
			return self.Lax_Wendroff_diff_Dx(si_data, coef)
		else :
			print("The func type name '{0}' of Dx is not defined.".format(sol_name))
			

class SiData():
	#data=np.array([])
	#old_data=np.array([])
	
	def set_data(self, delta_t, delta_x, data):
		self.dt = delta_t
		self.dx = delta_x
		self.N = len(data)
		self.data = np.array([])
		self.time = 0.0
		for i in range(0, self.N):
			self.data=np.append(self.data, data[i])
	
	def Dt(self, sol_name):
		return CalDiff().Dt(sol_name, self)
	
	#argsは指定された物理量にかかる係数:c
	def Dx(self, sol_name):
		return CalDiff().Dx(sol_name, self)
	
	def boundary(self):
		n=self.N
		self.data[0]=self.data_old[0]
		self.data[n-1]=self.data_old[n-1]
		
	def update(self, data, time):
		self.data_old = self.data
		self.data = data
		self.time = time

	def test(self, *args):
		for element in args :
			print("test: {0}\n".format(element))
		print(args)
		print(args[0])
		test_data = data
		print(test_data)


data=[5,6,7]
q=SiData()
#q.test(1)
#q.test(2,3,4)
#q.test(1,2,3,4,data)

#初期化
#data=[1,1,0,0]
#dx = 0.1
#dt = 1
#time = dt
#q=SiData()
#q.set_data(dt, dx, data)
#print("time : ",q.time,"data:",q.data)
#計算
#qn = (q.Dt("forward_diff") + q.Dx("centor_diff")).solve() #Dxの係数を変更したい
#q.update(qn,time)
#q.boundary()
#結果表示
#print("time : ",q.time,"data:",q.data)
