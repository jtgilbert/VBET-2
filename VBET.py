from PyQt5 import QtGui
from PyQt5.QtWidgets import QMainWindow, QFileDialog, QApplication, QDialog
import sys
import vbet_ui
import run_VBET


class VBETTool(QDialog, vbet_ui.Ui_Dialog):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)

        self.pushButton_stream.clicked.connect(lambda: self.file_browser(self.lineEdit_stream))
        self.pushButton_DEM.clicked.connect(lambda: self.file_browser(self.lineEdit_DEM))
        self.pushButton_scratch.clicked.connect(lambda: self.folder_browser(self.lineEdit_scratch))
        self.pushButton_da.clicked.connect(lambda: self.file_browser(self.lineEdit_da))
        self.pushButton_output.clicked.connect(lambda: self.file_save(self.lineEdit_output))

        self.buttonBox.accepted.connect(self.vbet)
        #self.buttonBox.rejected.connect(self.reject)

    def file_browser(self, txtControl):
        filename = QFileDialog.getOpenFileName(self, 'Open the File', '', 'Shapefiles (*.shp);; Rasters (*.tif);;',
                                               None, QFileDialog.DontUseNativeDialog)
        txtControl.setText(filename[0])

    def folder_browser(self, txtControl):
        folderpath = str(QFileDialog.getExistingDirectory(self, 'Select Folder', '', QFileDialog.DontUseNativeDialog))
        txtControl.setText(folderpath)

    def file_save(self, txtControl):
        filename = QFileDialog.getSaveFileName(self, 'Output File', '', 'Shapefile (*.shp)', None, QFileDialog.DontUseNativeDialog)
        txtControl.setText(filename[0])

    def vbet(self):
        inst = run_VBET.RunVBET()
        inst.params['network'] = str(self.lineEdit_stream.text())
        inst.params['dem'] = str(self.lineEdit_DEM.text())
        inst.params['out'] = str(self.lineEdit_output.text())
        inst.params['scratch'] = str(self.lineEdit_scratch.text())
        inst.params['lg_da'] = float(self.lineEdit_lgda.text())
        inst.params['med_da'] = float(self.lineEdit_medda.text())
        inst.params['lg_slope'] = self.SpinBox_lgslope.value()
        inst.params['med_slope'] = self.SpinBox_medslope.value()
        inst.params['sm_slope'] = self.SpinBox_smslope.value()
        inst.params['lg_buf'] = float(self.lineEdit_lgbuf.text())
        inst.params['med_buf'] = float(self.lineEdit_medbuf.text())
        inst.params['sm_buf'] = float(self.lineEdit_smbuf.text())
        inst.params['min_buf'] = float(self.lineEdit_minbuf.text())
        inst.params['dr_area'] = str(self.lineEdit_da.text())
        if str(self.lineEdit_exda.text()) != '':
            inst.params['da_field'] = str(self.lineEdit_exda.text())
        else:
            inst.params['da_field'] = None
        inst.params['lg_depth'] = self.SpinBox_lgdepth.value()
        inst.params['med_depth'] = self.SpinBox_meddepth.value()
        inst.params['sm_depth'] = self.SpinBox_smdepth.value()

        inst.run()


def main():
    app = QApplication(sys.argv)
    form = VBETTool()
    form.show()
    app.exec_()


if __name__ == '__main__':
    main()
