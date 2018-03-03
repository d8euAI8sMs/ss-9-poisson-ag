// poisson-agDlg.h : header file
//

#pragma once

#include <util/common/gui/SimulationDialog.h>
#include <util/common/gui/PlotControl.h>

#include "model.h"

// CPoissonAGDlg dialog
class CPoissonAGDlg : public CSimulationDialog
{
// Construction
public:
    CPoissonAGDlg(CWnd* pParent = NULL);    // standard constructor

// Dialog Data
    enum { IDD = IDD_POISSONAG_DIALOG };

    protected:
    virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support


// Implementation
protected:
    HICON m_hIcon;
    model::model_data data;

    // Generated message map functions
    virtual BOOL OnInitDialog();
    afx_msg void OnPaint();
    afx_msg HCURSOR OnQueryDragIcon();
    DECLARE_MESSAGE_MAP()
public:
    CPlotControl my_plot;
    void OnSimulation();
    afx_msg void OnBnClickedButton1();
    afx_msg void OnBnClickedButton2();
    afx_msg void OnBnClickedCheck3();
    BOOL m_bPointsVisible;
    BOOL m_bTriangulationVisible;
    BOOL m_bDirichletCellsVisible;
    BOOL m_bIsolinesVisible;
    BOOL m_bFieldLinesVisible;
};
