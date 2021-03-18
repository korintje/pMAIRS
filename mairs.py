import numpy as np
import pandas as pd

"""
angle: Angles of the substrate (experimental value) / degree.
theta: Tilt angle of the transition dipolar memoment calculated by OP and IP / degree.
"""

class MultiBeamsSet():

    def __init__(self):
        self.multibeams_bkg = MultiBeams()
        self.multibeams_smp = MultiBeams()
    
    def load_sample(self, beams):
        self.multibeams_smp = beams
    
    def load_background(self, beams):
        self.multibeams_bkg = beams

    def get_op_ip(self, ppolar=True):
        MAIRS = pd.DataFrame(columns=["wavenumber", "IP", "OP"])
        mairs_bkg = self.multibeams_bkg.get_op_ip(ppolar=ppolar)
        mairs_smp = self.multibeams_smp.get_op_ip(ppolar=ppolar)
        MAIRS["IP"] = -np.log(mairs_smp["IP"] / mairs_bkg["IP"])
        MAIRS["OP"] = -np.log(mairs_smp["OP"] / mairs_bkg["OP"])
        MAIRS["wavenumber"] = mairs_smp["wavenumber"][:]
        return MAIRS
    
    def get_thetas(self, n4H=1):
        mairs = self.get_op_ip()
        thetas = pd.DataFrame({
            "wavenumber": mairs["wavenumber"],
            "theta": np.rad2deg(np.arctan(np.sqrt(2 * mairs["IP"] / mairs["OP"] / n4H)))
        }) 
        return thetas


class MultiBeams():

    def __init__(self):
        self.beams = []

    def load_beam(self, angle, spectrum):
        wavenumbers = spectrum[0]
        transmittances = spectrum[1]
        self.beams.append(SingleBeam(angle, wavenumbers, transmittances))

    def get_op_ip(self, ppolar=True):
        """ By default, p-polarization is set True (pMAIRS) """
        MAIRS = pd.DataFrame(columns=["wavenumber", "IP", "OP"])
        beam_num = len(self.beams)
        if beam_num < 3:
            raise Exception(
                    "Error. At least 3 single-beams are required to calc OP/IP."
                    )

        # R行列（J. Phys. Chem. B 2002, 106, 4112.）
        _R = []
        f = 0 if ppolar is True else 1
        const = (4.0 / np.pi)**2.0
        for beam in self.beams:
            angle = np.deg2rad(beam.angle)
            vec = np.array([f + np.cos(angle)**2 + (np.sin(angle)**2 * np.tan(angle)**2), np.tan(angle)**2])
            _R.append(vec)
        R = const * np.array(_R)
        
        # 転置行列R^T
        R_trans = R.transpose()

        # 測定ビームスペクトルS
        S = np.array([beam.transmittances for beam in self.beams])

        # Compromise solution: 多次元空間における最小二乗法
        s_calc = np.linalg.inv(R_trans.dot(R)).dot(R_trans).dot(S)

        # Results
        MAIRS["IP"] = s_calc[0]
        MAIRS["OP"] = s_calc[1]
        MAIRS["wavenumber"] = self.beams[0].wavenumbers
        return MAIRS


class SingleBeam():

    def __init__(self, angle, wavenumbers, transmittances):
        self.angle = angle
        self.wavenumbers = wavenumbers
        self.transmittances = transmittances
