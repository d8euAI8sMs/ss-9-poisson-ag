// poisson-agDlg.h : header file
//

#pragma once

#include <util/common/gui/SimulationDialog.h>
#include <util/common/gui/PlotControl.h>

#include "model.h"
#include "afxwin.h"

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
    model::parameters params;

    enum _simmode {
        simmode_default,
        simmode_rev,
        simmode_dipole1_calc,
        simmode_dipole2_calc
    } simmode;

    struct _udata
    {
        CEdit U1, U2, Ud;
        CEdit d1x, d1y, d2x, d2y;
        geom::point2d_t d1, d2;
    } udata;

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
    UINT m_nIsolineCount;
    double m_fpIsolineDelta;
    UINT m_nFieldLineCount;
    geom::point2d_t DipoleAt(geom::point2d_t p, geom::point2d_t o);
    CEdit m_dip;
    afx_msg void OnCbnSelchangeCombo1();
    CComboBox m_simmode;
};
