FUNCTION rdtab, fname
;+
; NAME:
;       rdtab
; PURPOSE:
;       Reads a formatted table into a two dimensional array.
; CALLING SEQUENCE:
;       table = rdtab(fname)
; INPUTS:
;       fname:  Name of input text file.
; OUTPUTS:
;       table: A two dimensional array.
; MODIFICATION HISTORY:
;       BBT, September, 1994
;-
on_error, 2
on_ioerror, io_err
;
IF (n_elements(fname) eq 0) THEN message, 'Usage: table = rdtab(fname)'
;
;
; Open input file, and check for errors
openr, unit, fname, error=err, /get_lun
IF (err ne 0) THEN goto, finish
;
; Read the first record of data as text
text = ''
readf, unit, text
str = strcompress(strtrim(text, 2))
IF (strlen(str) eq 0) THEN message, 'The first record is empty'
;
; Determine the number of columns in the first record
ncols = 1
pos = strpos(str, ' ')
WHILE (pos ge 0) DO BEGIN
  pos = pos + 1
  ncols = ncols + 1
  pos = strpos(str, ' ', pos)
ENDWHILE
print, ncols
;
row   = fltarr(ncols)
table = row
;
; Read the first data record from string
reads, str, table
;
; Read and append the remaining rows to table
WHILE not eof(unit) DO BEGIN
  readf, unit, row
  table = [[table], [row]]
ENDWHILE
close, unit
;
goto, finish
io_err:
message, !err_string, /noname, /ioerror
finish:
on_ioerror, null
free_lun, unit
return, table
END
