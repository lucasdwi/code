function changeFont($dir)
# Goes through directory of .eps files and replaces all instances of "/Helvetica " with "/Arial " in order to maintain correct font for images exported from Matlab.

# INPUTS
# $dir = directory to search through
{
	cd $dir 
    $files = get-childitem -filter *.eps
	for ( $i = 0; $i -lt $files.Count; $i++) {
		(get-content $files[$i]) | foreach-object { $_ -replace "/Helvetica ","/Arial "} | set-content $files[$i]
	}
}