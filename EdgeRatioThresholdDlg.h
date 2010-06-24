/*
 * The information in this file is
 * Copyright(c) 2007 Ball Aerospace & Technologies Corporation
 * and is subject to the terms and conditions of the
 * GNU Lesser General Public License Version 2.1
 * The license text is available from   
 * http://www.gnu.org/licenses/lgpl.html
 */

#ifndef EDGE_RATIO_THRESHOLD_DLG_H
#define EDGE_RATIO_THRESHOLD_DLG_H

#include <QtGui/QDialog>


class QDoubleSpinBox;

class EdgeRatioThresholdDlg : public QDialog
{
   Q_OBJECT

public:
   EdgeRatioThresholdDlg(QWidget* pParent, double t1, double t2, double t3); 


private slots:
   void setSmallThreshold(double t);
   void setMeidanThreshold(double t);
   void setLargeThreshold(double t);

public:
   QDoubleSpinBox* mSmallThresholdBox;
   QDoubleSpinBox* mMedianThresholdBox;
   QDoubleSpinBox* mLargeThresholdBox;
   double getSmallThreshold();
   double getMedianThreshold();
   double getLargeThreshold();

private:
   double mSmallThreshold;
   double mMedianThreshold;
   double mLargeThreshold;
};

#endif
