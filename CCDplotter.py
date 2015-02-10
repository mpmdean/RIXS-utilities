import wx, os, glob
import pdb
import numpy as np
import matplotlib
# comment out the following to use wx rather than wxagg
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

from processCCD_image import *

# Add display box showing header of file?

#correct zoom using rubbish/rectange_select.py code
#self.axes.remove_patch(self.rect)
#somehow remove patch?

print "need to fix bug with sliding into end of range"

class ControlFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, title='File Browser', size=(400, 500), pos=(50, 50))

        panel = wx.Panel(self, -1)
        self.panel = panel
        bdr_sz = 10 # was 20
        
        # path TextCtrl
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.AddSpacer(bdr_sz)
        path_cap = wx.StaticText( panel, label='Search Term' )
        vbox.Add(path_cap, 0, wx.LEFT | wx.RIGHT | wx.ALIGN_LEFT, border=bdr_sz)
        self.path = wx.TextCtrl(panel, size=(140, 20), value = 
                      '/Users/mdean/Dropbox/ALS_feb_2015/data/2015 02 06/CCD Scan 2586/*.fits',
                      style=wx.TE_PROCESS_ENTER)
        self.path.Bind(wx.EVT_TEXT_ENTER, self.updatepath)
        vbox.Add(self.path, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, border=bdr_sz)
        vbox.AddSpacer(bdr_sz)
        
        # LIST BOX SHOWING FILES
        self.path_list, self.file_list = self.get_files(self.path.GetValue())
        listbox_cap = wx.StaticText( panel, label='Directory contents' )
        vbox.Add(listbox_cap, 0, wx.LEFT | wx.RIGHT | wx.ALIGN_LEFT, border=bdr_sz)
        self.listbox = wx.ListBox(panel, -1, wx.DefaultPosition, (170, 130), self.file_list, wx.LB_SINGLE)
        self.listbox.Bind(wx.EVT_LISTBOX, self.plotit)   # -- new binding!
        
        # buttons        
        vbox.Add(self.listbox, 1, wx.EXPAND | wx.LEFT |wx.RIGHT | wx.BOTTOM, border=bdr_sz)        
        buttonsize=(80, 20)
        hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.autoscale_but = wx.Button(panel, -1, 'Autoscale', size=buttonsize)
        self.Bind(wx.EVT_BUTTON, self.autoscaleit, self.autoscale_but)
        hbox.Add(self.autoscale_but, 0, wx.ALIGN_LEFT, border=bdr_sz)

#        self.vmin_TextCtrl = wx.TextCtrl(panel, -1, 'Autoscale', size=buttonsize)
#        self.Bind(wx.EVT_BUTTON, self.vminsetit, self.vmin_TextCtrl)
#        hbox.Add(self.vmin_TextCtrl, 0, wx.ALIGN_LEFT, border=bdr_sz)
        
        # sliders
        vbox.Add(hbox, 0, wx.LEFT | wx.RIGHT | wx.ALIGN_LEFT, border=bdr_sz)
        vbox.AddSpacer(bdr_sz)
        hbox = wx.BoxSizer(wx.HORIZONTAL) # added!
        self.vminslider = wx.Slider(panel, value=250, minValue=0, maxValue=500, 
                                    style=wx.SL_TOP | wx.SL_AUTOTICKS | wx.SL_LABELS)
        self.vminslider.SetTickFreq(50)
        self.Bind(wx.EVT_SCROLL, self.vminslideit, self.vminslider)
        hbox.Add(self.vminslider, 0, wx.LEFT | wx.RIGHT | wx.ALIGN_LEFT, border=bdr_sz)
        hbox.AddSpacer(bdr_sz)  # new
        
        self.vmaxslider = wx.Slider(panel, value=250, minValue=0, maxValue=500, 
                                    style=wx.SL_TOP | wx.SL_AUTOTICKS | wx.SL_LABELS)
        self.vmaxslider.SetTickFreq(50)
        self.Bind(wx.EVT_SCROLL, self.vmaxslideit, self.vmaxslider)
        hbox.Add(self.vmaxslider, 0, wx.LEFT | wx.RIGHT | wx.ALIGN_LEFT, border=bdr_sz)

        vbox.Add(hbox, 0, wx.LEFT | wx.RIGHT | wx.ALIGN_LEFT, border=bdr_sz) # new
        vbox.AddSpacer(bdr_sz)
        
        # Initiate everything
        panel.SetSizer(vbox)
        #self.Centre()
        self.Show(True)

        # Call plotframe        
        self.plotframe = PlotFrame(self)
        self.plotframe.Show()
        self.vmin = None
        self.vmax = None
        
    
    # functions to bind to events
    
    def updatepath(self, evt):
        print "update path"
        self.path_list, self.file_list = self.get_files(self.path.GetValue())
        self.listbox.Clear()
        for filelist in self.file_list:
            self.listbox.Append(filelist)
    
    def loadit(self):
        filepath = self.path_list[self.listbox.GetSelection()]
        self.CCD = CCD(fname_list=[filepath])

    def plotit(self, evt):
        self.loadit()
        #pdb.set_trace()
        self.M = self.CCD.raw_images[0]
        if (self.vmin == None) or (self.vmax == None):
            self.vmin = np.percentile(self.M,1)
            self.vmax = np.percentile(self.M,99)
            vrange = np.abs(self.vmax - self.vmin)
            self.vminslider.SetMin(self.vmin - 0.2 * vrange)
            self.vmaxslider.SetMin(self.vmin - 0.2 * vrange)
            self.vminslider.SetMax(self.vmax + 0.2 * vrange)
            self.vmaxslider.SetMax(self.vmax + 0.2 * vrange)
        
        self.vminslider.SetValue(self.vmin)
        self.vmaxslider.SetValue(self.vmax)
        
        self.plotframe.figure.clf()
        self.axes = self.plotframe.figure.add_subplot(111) 
        img = self.plotframe.figure.axes[0].imshow(self.M, vmin=self.vmin, vmax=self.vmax, aspect='auto', interpolation='none')            
        self.colorbarobj = plt.colorbar(img, ax=self.plotframe.figure.axes[0])
        self.plotframe.figure.canvas.draw()

    def autoscaleit(self, evt):
        print "Autoscale pressed"
        try:
            self.vmin = None
            self.vmax = None
            self.plotit(None)
        except AttributeError:
            "No image available"
            pass
    
    def vminslideit(self, evt):
        print "vmin slide"
        self.vmin = self.vminslider.GetValue()
        if self.vmin == self.vminslider.GetMin() or self.vmin == self.vminslider.GetMax():
            print "request number"
            dlg = wx.TextEntryDialog(self.panel, 'Enter value for minimum on colorbar',
                                     "Colorbar Dialog", "",  style=wx.OK)
            dlg.ShowModal()
            #pdb.set_trace()
            print dlg.GetValue()
            if dlg.GetValue().isnumeric():
                self.vmin = float(dlg.GetValue())
                self.vminslider.SetMin(self.vmin - np.abs(self.xmin * 0.1) - 1)
                self.vminslider.SetMax(self.vmin + - np.abs(self.vmin * 0.1) + 1)        
            dlg.Destroy()
        self.plotit(None)

    def vminsetit(self, evt):
        print "vmin set"
#        self.vmin = self.vminslider.GetValue()
#        self.plotit(None)

    def vmaxslideit(self, evt):
        print "vmax slide"
        self.vmax = self.vmaxslider.GetValue()
        if self.vmax == self.vmaxslider.GetMin() or self.vmax == self.vmaxslider.GetMax():
            #pdb.set_trace()
            print "request number"
            dlg = wx.TextEntryDialog(self.panel, 'Enter value for maximum on colorbar',
                                     "Colorbar Dialog", "",  style=wx.OK)
            dlg.ShowModal()
            print dlg.GetValue()
            if dlg.GetValue().isnumeric():
                self.vmax = float(dlg.GetValue())
                self.vmaxslider.SetMin(self.vmax - np.abs(self.vmax * 0.1) - 1)
                self.vmaxslider.SetMax(self.vmax + - np.abs(self.vmax * 0.1) + 1)            
            dlg.Destroy()
        self.plotit(None)
       

    def get_files(self, searchterm):
        path_list = glob.glob(searchterm)
        
        file_list = []
        for path in path_list:
            file_list.append(os.path.split(path)[1])
            
        return path_list, file_list
    
        

class PlotFrame(wx.Frame):
    def __init__(self, parent):
        wx.Frame.__init__(self, None, size=(300,500), pos=(500, 50), title='Plot Frame')
        self.parent = parent
        
        # initialize plot        
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)         
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        
        # Connect the mouse events to their relevant callbacks
        self.canvas.mpl_connect('button_press_event', self._onPress)
        self.canvas.mpl_connect('button_release_event', self._onRelease)
        self.canvas.mpl_connect('motion_notify_event', self._onMotion)
        
        self.pressed = False
        
        # Initialise the rectangle
        self.rect = Rectangle((0,0), 1, 1, facecolor='None', visible=False,
                                  edgecolor='k', linestyle='dashed')
        self.x0 = 0#None
        self.y0 = 0#None
        self.x1 = 0#None
        self.y1 = 0#None
        self.axes.add_patch(self.rect)
        print "shoudl add to axes[0]?"
        self.Fit()
        
    def _onPress(self, event):
        ''' Callback to handle the mouse being clicked and held over the canvas'''
        print "onPress"
        # Check the mouse press was actually on the canvas
        if event.xdata is not None and event.ydata is not None:
            print "In press event data"
            # Upon initial press of the mouse record the origin and record the mouse as pressed
            self.pressed = True
            self.rect.set_visible = True
            self.x0 = event.xdata
            self.y0 = event.ydata
            
            # Set the width and height and draw the rectangle
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
            #self.axes.add_patch(self.rect)
            self.canvas.draw()


    def _onRelease(self, event):
        '''Callback to handle the mouse being released over the canvas'''
        print "onRelease"
        # Check that the mouse was actually pressed on the canvas to begin with and this isn't a rouge mouse 
        # release event that started somewhere else
        if self.pressed:

            # Upon release draw the rectangle as a solid rectangle
            self.pressed = False

            # Check the mouse was released on the canvas, and if it wasn't then just leave the width and 
            # height as the last values set by the motion event
            if event.xdata is not None and event.ydata is not None:
                self.x1 = event.xdata
                self.y1 = event.ydata

            # Set the width and height and origin of the bounding rectangle
            self.boundingRectWidth =  self.x1 - self.x0
            self.boundingRectHeight =  self.y1 - self.y0
            self.bouningRectOrigin = (self.x0, self.y0)

            # Draw the bounding rectangle
            self.rect.set_width(self.boundingRectWidth)
            self.rect.set_height(self.boundingRectHeight)
            self.rect.set_xy((self.x0, self.y0))
            
            # set zoom
            xmin = min([self.x0, self.x1])
            xmax = max([self.x0, self.x1])
            ymin = min([self.y0, self.y1])
            ymax = max([self.y0, self.y1])
            self.figure.axes[0].axis([xmin, xmax, ymin, ymax])
            self.rect.set_visible = False

            self.canvas.draw()


    def _onMotion(self, event):
        '''Callback to handle the motion event created by the mouse moving over the canvas'''
        #print "onMotion"
        # If the mouse has been pressed draw an updated rectangle when the mouse is moved so 
        # the user can see what the current selection is
        if self.pressed:

            # Check the mouse was released on the canvas, and if it wasn't then just leave the width and 
            # height as the last values set by the motion event
            if event.xdata is not None and event.ydata is not None:
                self.x1 = event.xdata
                self.y1 = event.ydata
            
            # Set the width and height and draw the rectangle
            self.rect.set_visible = True
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
            #self.rectpatch = self.axes.add_patch(self.rect)
            self.canvas.draw()
        
    def OnPaint(self, event):
        self.canvas.draw()

if __name__ == "__main__":
    App=wx.PySimpleApp()
    ControlFrame().Show()
    App.MainLoop()