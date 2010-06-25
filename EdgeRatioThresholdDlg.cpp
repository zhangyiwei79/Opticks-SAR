/*
 * The information in this file is
 * Copyright(c) 2007 Ball Aerospace & Technologies Corporation
 * and is subject to the terms and conditions of the
 * GNU Lesser General Public License Version 2.1
 * The license text is available from   
 * http://www.gnu.org/licenses/lgpl.html
 */

#include "AppAssert.h"
#include "AppVerify.h"
#include "EdgeRatioThresholdDlg.h"

#include <QtGui/QDoubleSpinBox>
#include <QtGui/QLabel>
#include <QtGui/QLayout>
#include <QtGui/QPushButton>

using namespace std;

EdgeRatioThresholdDlg::EdgeRatioThresholdDlg(QWidget* pParent, double smallThreshold, double medianThreshold, double largeThreshold) : QDialog(pParent),
   mSmallThresholdBox(NULL), mMedianThresholdBox(NULL),mLargeThresholdBox(NULL)
{
   setWindowTitle("Ratio Threshold");

   QGridLayout* pLayout = new QGridLayout(this);
   pLayout->setMargin(10);
   pLayout->setSpacing(5);

   QLabel* pLable1 = new QLabel("3 by 3 window threshold: ", this);
   pLayout->addWidget(pLable1, 0, 0);

   mSmallThresholdBox = new QDoubleSpinBox(this);
   mSmallThresholdBox->setRange(0.1, 1);
   mSmallThresholdBox->setSingleStep(0.01);
   mSmallThresholdBox->setValue(smallThreshold);
   mSmallThreshold = smallThreshold;
   pLayout->addWidget(mSmallThresholdBox, 0, 1, 1, 2);
   //pLayout->setColumnStretch(1, 10);

   QLabel* pLable2 = new QLabel("5 by 5 window threshold: ", this);
   pLayout->addWidget(pLable2, 1, 0);

   mMedianThresholdBox = new QDoubleSpinBox(this);
   mMedianThresholdBox->setRange(0.1, 1);
   mMedianThresholdBox->setSingleStep(0.01);
   mMedianThresholdBox->setValue(medianThreshold);
   mMedianThreshold = medianThreshold;
   pLayout->addWidget(mMedianThresholdBox, 1, 1, 1, 2);

   QLabel* pLable3 = new QLabel("7 by 7 window threshold: ", this);
   pLayout->addWidget(pLable3, 2, 0);

   mLargeThresholdBox = new QDoubleSpinBox(this);
   mLargeThresholdBox->setRange(0.1, 1);
   mLargeThresholdBox->setSingleStep(0.01);
   mLargeThresholdBox->setValue(largeThreshold);
   mLargeThreshold = largeThreshold;
   pLayout->addWidget(mLargeThresholdBox, 2, 1, 1, 2);

   //pLayout->setRowStretch(2, 10);

   QHBoxLayout* pRespLayout = new QHBoxLayout;
   pLayout->addLayout(pRespLayout, 3, 0, 1, 3);

   QPushButton* pAccept = new QPushButton("OK", this);
   pRespLayout->addStretch();
   pRespLayout->addWidget(pAccept);

   QPushButton* pReject = new QPushButton("Cancel", this);
   pRespLayout->addWidget(pReject);

   connect(pAccept, SIGNAL(clicked()), this, SLOT(accept()));
   connect(pReject, SIGNAL(clicked()), this, SLOT(reject()));
   
   connect(mSmallThresholdBox, SIGNAL(valueChanged(double)), this, SLOT(setSmallThreshold(double)));
   connect(mMedianThresholdBox, SIGNAL(valueChanged(double)), this, SLOT(setMedianThreshold(double)));
   connect(mLargeThresholdBox, SIGNAL(valueChanged(double)), this, SLOT(setLargeThreshold(double)));
}

void EdgeRatioThresholdDlg::setSmallThreshold(double t)
{
	mSmallThreshold = t;
}

void EdgeRatioThresholdDlg::setMeidanThreshold(double t)
{
	mMedianThreshold = t;
}

void EdgeRatioThresholdDlg::setLargeThreshold(double t)
{
	mLargeThreshold = t;
}

double EdgeRatioThresholdDlg::getSmallThreshold()
{
	return mSmallThreshold;
}

double EdgeRatioThresholdDlg::getMedianThreshold()
{
	return mMedianThreshold;
}

double EdgeRatioThresholdDlg::getLargeThreshold()
{
	return mLargeThreshold;
}



