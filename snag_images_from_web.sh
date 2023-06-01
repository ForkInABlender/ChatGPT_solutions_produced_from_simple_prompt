echo $1| while read -r url; do     curl -s "$url" | grep -o '<img[^>]*src="[^"]*"' | grep -o 'src="[^"]*"' | grep -o '"[^"]*"' | tr -d '"'; done
