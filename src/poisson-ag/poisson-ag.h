// poisson-ag.h : main header file for the PROJECT_NAME application
//

#pragma once

#ifndef __AFXWIN_H__
    #error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"        // main symbols


// CPoissonAGApp:
// See poisson-ag.cpp for the implementation of this class
//

class CPoissonAGApp : public CWinApp
{
public:
    CPoissonAGApp();

// Overrides
public:
    virtual BOOL InitInstance();

// Implementation

    DECLARE_MESSAGE_MAP()
};

extern CPoissonAGApp theApp;
