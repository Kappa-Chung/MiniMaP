from tkinter import StringVar
from tkinter.filedialog import askopenfilenames, asksaveasfilename
import tkinter as tk
from tkinter import ttk
import numpy as np
import pandas as pd
import netCDF4 as nc
import math
import pickle
import os
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.ticker import FuncFormatter
from pathlib import Path
import matplotlib.pyplot as plt
from PIL import Image, ImageTk

class Processing():
    def __init__(self):
        pass
    def open_file(self):
        filename = askopenfilenames(title='Select a file',filetypes=[('CDF files','*.cdf')])
        self.filename = filename
        a = []
        for file in self.filename:
            a.append(file.split('/')[-1][:-4])
        self.filename_list = a
        a = pd.DataFrame(a)
        path.set(a.to_string(index = None, header = None))
    def machine_learning(self):
        spectrum_final_output = self.spectra_preprocessing()
        self.spectrum = spectrum_final_output
        biopsy_prediciton = self.prediction(spectrum_final_output)
        self.label_output(biopsy_prediciton) 
    def correction(self,mass1,mass2,int_array_2):
        index1 = np.where(np.isclose(int_array_2[0,:], mass1))
        index2 = np.where(np.isclose(int_array_2[0,:], mass2))
        correct_mass = max(int_array_2[1, index1[0][0]:index2[0][0]])
        index = np.where(np.isclose(int_array_2[1,:], correct_mass))
        index = index[0][0]
        return index

    def spectra_preprocessing(self,bin_size=1,z=1): 
        bin_parameter = int(bin_size/0.05)
        for file in self.filename:
            data = nc.Dataset(file,'a')
            inten = data.variables['intensity_values'][:]
            mass = data.variables['mass_values'][:]
            TIC = sum(data.variables['total_intensity'][:])
            ms_list = np.arange(np.min(mass), np.max(mass), 0.05)
            ms_list = np.around(ms_list,2)
            mass = np.around(mass,2)
            
            spectrum = np.zeros((2,len(ms_list)))
            spectrum[0,:] = ms_list
            for i in range(len(ms_list)):
                ind = np.where(mass == ms_list[i])
                spectrum[1, i] = sum(inten[ind])
                
            ref1 = self.correction(518,519,spectrum)
            ref2 = self.correction(725,726.5,spectrum)
            ref3 = self.correction(782,784,spectrum)
            ref4 = self.correction(808,809.5,spectrum)
            
            Obs_peak=[ref1,ref2,ref3,ref4]
            Cal_peak=[518.32,725.56,782.57,808.58]
            
            # Calculate average mass error
            error=(sum(ms_list[Obs_peak])-sum(Cal_peak))/4
            # Calculate number of position to move
            error_parameter=np.ceil(error/0.05)

            if error_parameter > 0:
                # m/z is higher than expected.
                spectrum_corrected = np.zeros((2,spectrum.shape[1]))
                spectrum_corrected[0,:-int(error_parameter)] = ms_list[:-int(error_parameter)]
                spectrum_corrected[1,:-int(error_parameter)] = spectrum[1,int(error_parameter):]
    
            elif error_parameter < 0:
                # m/z is lower than expected.
                error_parameter = error_parameter*(-1)
                spectrum_corrected = np.zeros((2,spectrum.shape[1]+int(error_parameter)))
                spectrum_corrected[0,:-int(error_parameter)] = ms_list[:]
                spectrum_corrected[1,int(error_parameter):] = spectrum[1,:]
    
            elif error_parameter == 0:
                # No calibration required
                spectrum_corrected = spectrum
            
            ms_list_bin = np.arange(math.ceil(min(ms_list))+1.5,math.floor(max(ms_list))-3.5,bin_size)
            spectrum_final = np.zeros((2,len(ms_list_bin)))
            spectrum_final[0,:] = ms_list_bin    
            cal_range = int(bin_parameter/2)
            
            for i in range(len(ms_list_bin)):
                index = np.where(np.isclose(ms_list_bin[i],spectrum_corrected[0,:]))
                index = index[0][0]
                spectrum_final[1,i] = sum(spectrum_corrected[1,index-cal_range:index+cal_range])/TIC
               
            if z == 1:
                spectrum_final_output = spectrum_final
            if z != 1:
                spectrum_final_output = np.vstack([spectrum_final_output,spectrum_final[1,:]])
            z+=1
        return spectrum_final_output
        
    def prediction(self,spectrum_final_output):
        script_dir = Path('Dr_Hsu_MiniMaP.py').parent.absolute()
        svc = pickle.load(open(os.path.join(script_dir,'Model.pkl'), 'rb'))
        biopsy_prediction = svc.predict(spectrum_final_output[1:,:])
        return biopsy_prediction
    
    def label_output(self,biopsy_prediction):
        biopsy_prediction = pd.DataFrame(list(biopsy_prediction))
        biopsy_prediction = biopsy_prediction.replace(0,'Benign')
        biopsy_prediction = biopsy_prediction.replace(1,'Malignant')
        prediction_result.set(biopsy_prediction.to_string(index=False,header=False, justify='left'))

    def plot_spectra(self, idx, mass_start, mass_end):
        ms_list = self.spectrum[0,:]
        spectrum = self.spectrum[1+idx,:]
        spectrum -= min(spectrum)

        mass_start, mass_end = float(mass_start), float(mass_end)
        mass_start += 0.5
        mass_end -= 0.5
 
        fig, ax = plt.subplots(figsize=(9,3))
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
        a, b = np.where(ms_list == mass_start)[0][0], np.where(ms_list == mass_end)[0][0]
        plt.plot(ms_list[a:b], spectrum[a:b]/max(spectrum[a:b]), color='k')
        plt.xlim([mass_start-5, mass_end+5])
        plt.ylim([0, 1])
        plt.xlabel('m/z', style='italic', fontname='Times New Roman')
        plt.ylabel('Relative Abundance', fontname='Times New Roman')

        def to_percent(temp, position):
          return '%1.0f'%(100*temp) + '%'
        plt.gca().yaxis.set_major_formatter(FuncFormatter(to_percent))
        self.fig = fig
                 
    def photo_resize(self, path, ratio):
        photo = Image.open(path)
        x, y = photo.size
        photo = photo.resize((int(x*ratio), int(y*ratio)))
        photo = ImageTk.PhotoImage(photo)
        return photo
    
    def export_spectra(self):
        filepath = asksaveasfilename(confirmoverwrite = True, 
                                     title = 'Save file', 
                                     filetypes = [('png files','*.png')])
        plt.savefig(filepath)
        
        
class GUI():
    def __init__(self):   
        self.section_borderline_button_func()
        self.section_borderline_button_2_func()
        self.section_borderline_button_3_func()
        self.section_borderline_button_4_func()
        self.file_open_label_func()
        self.file_open_blank_func()
        self.prediction_result_label_func()
        self.prediction_result_blank_func()
        self.file_button_func()
        self.running_button_func()
        self.export_spectra_button_func()
        self.spectra_combobox_list_func()
        self.ms_profile_label_func()
        self.mass_range_label_func()
        self.mass_range_func()
        self.display_spectra_button_func()
        self.display_spectra_diagram_func()
        self.trademark_label_func()
        self.background_color = '#FFFFBB'
        
    def file_open_label_func(self):

        self.file_open_label = ttk.Label(root)
        self.file_open_label['text'] = 'Open File'
        self.file_open_label['font'] = ('Calibri', 16)
        self.file_open_label.place(relx = 0.21, rely = 0.045)

    def file_open_blank_func(self):
        self.file_open_blank = ttk.Label(root)
        self.file_open_blank['textvariable'] = path
        self.file_open_blank['width'] = 40
        self.file_open_blank['font'] = ('Calibri', 12)
        self.file_open_blank['anchor'] = 'center'
        self.file_open_blank['background'] = 'white'
        self.file_open_blank['borderwidth'] = 1
        self.file_open_blank['relief'] = 'solid'
        self.file_open_blank.place(relx = 0.21, rely = 0.13)
        
    def prediction_result_label_func(self):
        self.prediction_result_label = ttk.Label(root)
        self.prediction_result_label['text'] = 'Prediction Result'
        self.prediction_result_label['font'] = ('Calibri', 16)
        self.prediction_result_label.place(relx = 0.63, rely = 0.045)

    def prediction_result_blank_func(self):
        self.prediction_result_blank = ttk.Label(root)
        self.prediction_result_blank['width'] = 20
        self.prediction_result_blank['textvariable'] = prediction_result
        self.prediction_result_blank['font'] = ('Calibri', 12)
        self.prediction_result_blank['background'] = 'white'
        self.prediction_result_blank['anchor'] = 'center'
        self.prediction_result_blank['borderwidth'] = 1
        self.prediction_result_blank['relief'] = 'solid'
        self.prediction_result_blank.place(relx = 0.63, rely = 0.13)
        
    def file_button_func(self):
        self.file_button = ttk.Button(root)
        self.file_button['text'] = 'Open file'
        self.file_button['image'] = photo
        self.file_button['command'] = P.open_file
        self.file_button.place(relx = 0.035, rely = 0.03)
        
    def running_button_func(self):
        self.running_button = ttk.Button(root)
        self.running_button['text'] = 'Run'
        self.running_button['image'] = photo2
        self.running_button['command'] = P.machine_learning
        self.running_button.place(relx = 0.035, rely = 0.28)
    def export_spectra_button_func(self):
        self.export_spectra_button = ttk.Button(root)
        self.export_spectra_button['text'] = 'Export spectra'
        self.export_spectra_button['image'] = photo3
        self.export_spectra_button['command'] = P.export_spectra
        self.export_spectra_button.place(relx = 0.035, rely = 0.53)
        
    def mass_range_label_func(self):
        self.mass_range_label = ttk.Label(root)
        self.mass_range_label['text'] = 'Mass Range (m/z)'
        self.mass_range_label['font'] = ('Calibri', 16)
        self.mass_range_label.place(relx = 0.63, rely = 0.37)
        
    def ms_profile_label_func(self):
        self.mass_range_label = ttk.Label(root)
        self.mass_range_label['text'] = 'MS Profile'
        self.mass_range_label['font'] = ('Calibri', 16)
        self.mass_range_label.place(relx = 0.21, rely = 0.37)
        
    def mass_range_func(self):
        self.mass_range_start = ttk.Entry(root)
        self.mass_range_start['width'] = 10
        self.mass_range_start['textvariable'] = mass_start
        self.mass_range_start['justify'] = 'center'
        self.mass_range_start['font'] = ('Calibri', 10)
        self.mass_range_start.place(relx = 0.63, rely = 0.44)
        
        self.mass_range_end = ttk.Entry(root)
        self.mass_range_end['width'] = 10
        self.mass_range_end['textvariable'] = mass_end
        self.mass_range_end['justify'] = 'center'
        self.mass_range_end['font'] = ('Calibri', 10)
        self.mass_range_end.place(relx = 0.76, rely =0.44)
        
        self.mass_range_tilde = ttk.Label(root)
        self.mass_range_tilde['text'] = '～'
        self.mass_range_tilde['font'] = ('Calibri', 16)
        self.mass_range_tilde.place(relx = 0.72, rely =0.435)
    
    def spectra_combobox_list_func(self):
        self.spectra_combobox_list = ttk.Combobox(root)
        self.spectra_combobox_list['postcommand'] = self.spectra_combobox_list_value_func
        self.spectra_combobox_list['state'] = 'readonly'
        self.spectra_combobox_list['width'] = 40
        self.spectra_combobox_list['font'] = ('Calibri', 12)
        self.spectra_combobox_list['justify'] = 'center'
        self.spectra_combobox_list.place(relx = 0.21, rely = 0.435)

    def spectra_combobox_list_value_func(self):
        self.spectra_combobox_list['values'] = P.filename_list
        self.spectra_combobox_list['font'] = ('Calibri', 12)

    def display_spectra_button_func(self):
        self.display_spectra_button = tk.Button(root)
        self.display_spectra_button['text'] = 'Enter'
        self.display_spectra_button['width'] = 10
        self.display_spectra_button['height'] = 1
        self.display_spectra_button.bind('<ButtonRelease-1>', self.display_spectra_func)
        root.bind('<Return>', self.display_spectra_func)
        self.display_spectra_button.place(relx = 0.865, rely = 0.432)
        
    def display_spectra_func(self, event):    
        mass_start = self.mass_range_start.get()
        mass_end = self.mass_range_end.get()
        mass_start, mass_end = float(mass_start), float(mass_end)
        
        if mass_start < 501 or mass_end > 990:
            self.mass_range_error_messagebox = tk.messagebox.showerror('Error', 'mass range should be in 501 to 990' )
        else:
            a = P.filename_list
            b = self.spectra_combobox_list.get()
            idx = a.index(b)
            P.plot_spectra(idx, mass_start, mass_end)
            self.display_spectra = FigureCanvasTkAgg(P.fig, root)
            P.relx = 0.22
            P.rely = 0.51
            self.display_spectra.get_tk_widget().place(relx=P.relx, rely=P.rely)
     
    def display_spectra_diagram_func(self):
        self.display_spectra_diagram = tk.Label()
        self.display_spectra_diagram['width'] = 94
        self.display_spectra_diagram['background'] = 'white'
        self.display_spectra_diagram['height'] = 15
        self.display_spectra_diagram['borderwidth'] = 2
        self.display_spectra_diagram['relief'] = 'ridge'
        self.display_spectra_diagram.place(relx = 0.209, rely = 0.5)
    
    def section_borderline_button_func(self):
        self.section_borderline_button = tk.Label()
        self.section_borderline_button['width'] = 23
        self.section_borderline_button['height'] = 27
        self.section_borderline_button['borderwidth'] = 0
        self.section_borderline_button['relief'] = 'ridge'
        self.section_borderline_button.place(relx = 0.021, rely = 0.01)
        
    def section_borderline_button_2_func(self):
        self.section_borderline_button = tk.Label()
        self.section_borderline_button['width'] = 50
        self.section_borderline_button['height'] = 10
        self.section_borderline_button['borderwidth'] = 2
        self.section_borderline_button['relief'] = 'ridge'
        self.section_borderline_button.place(relx = 0.2, rely = 0.03)
        
    def section_borderline_button_3_func(self):
        self.section_borderline_button = tk.Label()
        self.section_borderline_button['width'] = 45
        self.section_borderline_button['height'] = 10
        self.section_borderline_button['borderwidth'] = 2
        self.section_borderline_button['relief'] = 'ridge'
        self.section_borderline_button.place(relx = 0.62, rely = 0.03)
        
    def section_borderline_button_4_func(self):
        self.section_borderline_button = tk.Label()
        self.section_borderline_button['width'] = 99
        self.section_borderline_button['height'] = 20
        self.section_borderline_button['borderwidth'] = 2
        self.section_borderline_button['relief'] = 'ridge'
        self.section_borderline_button.place(relx = 0.2, rely = 0.365)
        
    def trademark_label_func(self):
        self.trademark_label = tk.Label()
        self.trademark_label['text'] = 'Cheng-Chih, Hsu Lab \n Deparment of chemistry, NTU \n Tel +886-2-3366-1681'
        self.trademark_label['font'] = ('Calibri', 10)
        self.trademark_label.place(relx = 0.005, rely = 0.84)
        

if __name__ == '__main__':
    root = tk.Tk()
    root.title('Dr.Hsu')
    root.geometry('900x500')
    
    path = StringVar()
    prediction_result = StringVar()
    prediction_result_file = StringVar()
    mass_start = StringVar(root, value = 501)
    mass_end = StringVar(root, value = 990)
    script_dir = Path('Dr_Hsu_MiniMaP.py').parent.absolute()
    
    P = Processing()
    photo = P.photo_resize(script_dir/'icon'/'open_file.png', 0.2)
    photo2 = P.photo_resize(script_dir/'icon'/'prediction_2.png', 0.2)
    photo3 = P.photo_resize(script_dir/'icon'/'export_spectra.png', 0.2)
    G = GUI()

    tk.mainloop()