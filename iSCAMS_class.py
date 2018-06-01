# iSCAMS class file.
# Charles Hill
# 01/06/18

class iSCAMS:
    def __init__(self,Cf,c_range=(0.0,0.15),m_range=(0.0,1500),Mass=True,bin_type='knuth',order=3):
        self.contrast = abs(Cf)
        self.instances = len(Cf)
        self.mass_data = []
        self.popt = []
        self.p_guess = []
        self.x = []
        self.fit = []
        self.bins = bin_type
        self.Mass = Mass
        self.c_range = c_range
        self.m_range = m_range
        self.order = order

    def func(self, x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)

    def Fit_Gaussian(self):
        plt.figure()
        if self.Mass == True:
            n, bins, pathces = hist(self.mass_data,bins=self.bins,normed=True,align='mid')
        else:
            n, bins, pathces = hist(self.contrast,bins=self.bins,normed=True,align='mid')
        #print("Number of Bins:",len(bins))
        while True:
            try:
                self.popt, pcov = curve_fit(self.func, bins[:-1], n, p0=self.p_guess)
                break
            except RuntimeError:
                print("Error - curve_fit failed")
                self.p_guess = self.p_guess[:-3]

        
        if self.Mass == True:
            self.x = np.linspace(self.m_range[0],self.m_range[1],1500)
        else:
            self.x = np.linspace(self.c_range[0],self.c_range[1],1500)
        self.fit = self.func(self.x, self.*popt)

        plt.plot(self.x,self.fit, 'b--')
        plt.yticks([])
        if self.Mass == True:
            plt.xlabel("Mass (kDa)",fontsize = 12, color = 'blue')
        else:
            plt.xlabel("Contrast",fontsize = 12,color='blue')
        plt.ylabel("Probability Density",fontsize = 12,color='blue')
        plt.show()

    def Auto_Gauss(self):

        plt.figure
        if self.Mass == True:
            n, bins, pat = hist(self.mass_data,bins='knuth',align='left')
        else:
            n, bins, pat = hist(self.contrast,bins='knuth',align='left')
        plt.show()

        Rel_Max = argrelextrema(n,np.greater,order=3)
        ctr = bins[Rel_Max]
        amp = n[Rel_Max]
        Wid = np.zeros(len(n[Rel_Max]))
        j = 0
        for idx in Rel_Max[0]:
            i = 1
            while n[idx]/2.0 < n[idx+i]:
                i += 1
            ind = idx + i
            Wid[j] = (bins[ind]-bins[idx])*2
            j += 1

        p_guess = np.zeros(len(Wid)*3)
        j = 0
        for i in range(0,len(p_guess),3):
            p_guess[i] = ctr[j]
            p_guess[i+1] = amp[j]
            p_guess[i+2] = Wid[j]
            j += 1
        
        print("Auto p_guess:",p_guess)
        return p_guess

