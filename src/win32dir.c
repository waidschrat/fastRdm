/* win32dir.c - emulate opendir/readdir/closedir on MS Windows systems
 *
 * This is a modified version of the file build/win32/dirent/dirent.c
 * as distributed by glib-2.12.4.  The original file comes with the
 * following notice:
 *
 *   dirent.c
 *   This file has no copyright assigned and is placed in the Public
 *   Domain.  This file is a part of the mingw-runtime package.  No
 *   warranty is given; refer to the file DISCLAIMER within the
 *   package.
 *
 *   Derived from DIRLIB.C by Matt J. Weinstein
 *   This note appears in the DIRLIB.H
 *   DIRLIB.H by M. J. Weinstein   Released to public domain 1-Jan-89
 *
 *   Updated by Jeremy Bettis <jeremy@hksys.com>
 *   Significantly revised and rewinddir, seekdir and telldir added by
 *   Colin Peters <colin@fu.is.saga-u.ac.jp>
 *
 * All changed from the original file are copyright (C) 2006 Jochen
 * Voss and can be redistributed under the terms of the GNU General
 * Public License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 */

#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <io.h>
#include <direct.h>
#define WIN32_LEAN_AND_MEAN
#include <windows.h> /* for GetFileAttributes */
#include <tchar.h>

#include "win32dir.h"

#ifdef _UNICODE
#define _tdirent	_wdirent
#define _TDIR		_WDIR
#define _topendir	_wopendir
#define _tclosedir	_wclosedir
#define _treaddir	_wreaddir
#define _trewinddir	_wrewinddir
#define _ttelldir	_wtelldir
#define _tseekdir	_wseekdir
#else
#define _tdirent	dirent
#define _TDIR		DIR
#define _topendir	opendir
#define _tclosedir	closedir
#define _treaddir	readdir
#define _trewinddir	rewinddir
#define _ttelldir	telldir
#define _tseekdir	seekdir
#endif

#define SUFFIX	_T("*")
#define	SLASH	_T("\\")


/*
 * opendir
 *
 * Returns a pointer to a DIR structure appropriately filled in to begin
 * searching a directory.
 */
_TDIR *
_topendir (const _TCHAR *szPath)
{
	_TDIR *nd;
	unsigned int rc;
	_TCHAR szFullPath[MAX_PATH];

	errno = 0;

	if (!szPath)
		{
			errno = EFAULT;
			return (_TDIR *) 0;
		}

	if (szPath[0] == _T('\0'))
		{
			errno = ENOTDIR;
			return (_TDIR *) 0;
		}

	/* Attempt to determine if the given path really is a directory. */
	rc = GetFileAttributes (szPath);
	if (rc == (unsigned int)-1)
		{
			/* call GetLastError for more error info */
			errno = ENOENT;
			return (_TDIR *) 0;
		}
	if (!(rc & FILE_ATTRIBUTE_DIRECTORY))
		{
			/* Error, entry exists but not a directory. */
			errno = ENOTDIR;
			return (_TDIR *) 0;
		}

	/* Make an absolute pathname.  */
	_tfullpath (szFullPath, szPath, MAX_PATH);

	/* Allocate enough space to store DIR structure and the complete
	 * directory path given. */
	nd = (_TDIR *) malloc (sizeof (_TDIR) + (_tcslen(szFullPath) + _tcslen (SLASH) +
											 _tcslen(SUFFIX) + 1) * sizeof(_TCHAR));

	if (!nd)
		{
			/* Error, out of memory. */
			errno = ENOMEM;
			return (_TDIR *) 0;
		}

	/* Create the search expression. */
	_tcscpy (nd->dd_name, szFullPath);

	/* Add on a slash if the path does not end with one. */
	if (nd->dd_name[0] != _T('\0') &&
		nd->dd_name[_tcslen (nd->dd_name) - 1] != _T('/') &&
		nd->dd_name[_tcslen (nd->dd_name) - 1] != _T('\\'))
		{
			_tcscat (nd->dd_name, SLASH);
		}

	/* Add on the search pattern */
	_tcscat (nd->dd_name, SUFFIX);

	/* Initialize handle to -1 so that a premature closedir doesn't try
	 * to call _findclose on it. */
	nd->dd_handle = -1;

	/* Initialize the status. */
	nd->dd_stat = 0;

	/* Initialize the dirent structure. ino and reclen are invalid under
	 * Win32, and name simply points at the appropriate part of the
	 * findfirst_t structure. */
	nd->dd_dir.d_ino = 0;
	nd->dd_dir.d_reclen = 0;
	nd->dd_dir.d_namlen = 0;
	memset (nd->dd_dir.d_name, 0, FILENAME_MAX);

	return nd;
}


/*
 * readdir
 *
 * Return a pointer to a dirent structure filled with the information on the
 * next entry in the directory.
 */
struct _tdirent *
_treaddir (_TDIR * dirp)
{
	errno = 0;

	/* Check for valid DIR struct. */
	if (!dirp)
		{
			errno = EFAULT;
			return (struct _tdirent *) 0;
		}

	if (dirp->dd_stat < 0)
		{
			/* We have already returned all files in the directory
			 * (or the structure has an invalid dd_stat). */
			return (struct _tdirent *) 0;
		}
	else if (dirp->dd_stat == 0)
		{
			/* We haven't started the search yet. */
			/* Start the search */
			dirp->dd_handle = _tfindfirst (dirp->dd_name, &(dirp->dd_dta));

			if (dirp->dd_handle == -1)
				{
					/* Whoops! Seems there are no files in that
					 * directory. */
					dirp->dd_stat = -1;
				}
			else
				{
					dirp->dd_stat = 1;
				}
		}
	else
		{
			/* Get the next search entry. */
			if (_tfindnext (dirp->dd_handle, &(dirp->dd_dta)))
				{
					/* We are off the end or otherwise error.
					   _findnext sets errno to ENOENT if no more file
					   Undo this. */
					DWORD winerr = GetLastError();
					if (winerr == ERROR_NO_MORE_FILES)
						errno = 0;
					_findclose (dirp->dd_handle);
					dirp->dd_handle = -1;
					dirp->dd_stat = -1;
				}
			else
				{
					/* Update the status to indicate the correct
					 * number. */
					dirp->dd_stat++;
				}
		}

	if (dirp->dd_stat > 0)
		{
			/* Successfully got an entry. Everything about the file is
			 * already appropriately filled in except the length of the
			 * file name. */
			dirp->dd_dir.d_namlen = _tcslen (dirp->dd_dta.name);
			_tcscpy (dirp->dd_dir.d_name, dirp->dd_dta.name);
			return &dirp->dd_dir;
		}

	return (struct _tdirent *) 0;
}


/*
 * closedir
 *
 * Frees up resources allocated by opendir.
 */
int
_tclosedir (_TDIR * dirp)
{
	int rc;

	errno = 0;
	rc = 0;

	if (!dirp)
		{
			errno = EFAULT;
			return -1;
		}

	if (dirp->dd_handle != -1)
		{
			rc = _findclose (dirp->dd_handle);
		}

	/* Delete the dir structure. */
	free (dirp);

	return rc;
}
