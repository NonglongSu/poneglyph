
#Add .fa suffix to nameList. 
#awk 'NF{print $0 ".fa"}' inFile


while IFS= read -r f;do
	if ! [[ -e $2/$f ]]; then
		printf '%s is missing in %s\n' "$f" "$2"
	fi
done < "$1"

