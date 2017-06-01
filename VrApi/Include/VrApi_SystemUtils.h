/************************************************************************************

Filename    :   VrApi_SystemUtils.h
Content     :   Interface for SystemUtils functionality.
Created     :   August 15, 2014
Authors     :   Gloria Kennickell, Jonathan E. Wright

Copyright   :   Copyright 2014 Oculus VR, LLC. All Rights reserved.

*************************************************************************************/
#ifndef OVR_VrApi_SystemUtils_h
#define OVR_VrApi_SystemUtils_h

#include "VrApi_Config.h"
#include "VrApi_Types.h"

#if defined( __cplusplus )
extern "C" {
#endif

typedef enum
{
	VRAPI_SYS_UI_GLOBAL_MENU,					// Display the Universal Menu.
	VRAPI_SYS_UI_CONFIRM_QUIT_MENU,				// Display the 'Confirm Quit' Menu.

	VRAPI_SYS_UI_KEYBOARD_MENU,					// Display a Keyboard Menu for editing a single string.
	VRAPI_SYS_UI_FILE_DIALOG_MENU,				// Display a Folder Browser Menu for selecting the path to a file or folder.

} ovrSystemUIType;

// Display a specific System UI.
OVR_VRAPI_EXPORT bool vrapi_ShowSystemUI( const ovrJava * java, const ovrSystemUIType type );

// ----DEPRECATED
// This function is DEPRECATED. Please do not write any new code which
// relies on it's use.
// Display a specific System UI and pass extra JSON text data.
OVR_VRAPI_EXPORT bool vrapi_ShowSystemUIWithExtra( const ovrJava * java, const ovrSystemUIType type, const char * extraJsonText );

// Launch the Oculus Home application.
// NOTE: This exits the current application by executing a finishAffinity.
OVR_VRAPI_EXPORT void vrapi_ReturnToHome( const ovrJava * java );

// Display a Fatal Error Message using the System UI.
OVR_VRAPI_EXPORT void vrapi_ShowFatalError( const ovrJava * java, const char * title, const char * message,
		const char * fileName, const unsigned int lineNumber );

#if defined( __cplusplus )
}	// extern "C"
#endif

#endif	// OVR_VrApi_SystemUtils_h
