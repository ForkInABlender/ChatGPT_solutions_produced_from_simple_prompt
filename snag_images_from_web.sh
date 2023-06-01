# Dylan Kenneth Eliot+GPT+Dan-6.0-prompt

# Using a triple quoted string, read from each url only the image URLs from <img> tags

echo $1| while read -r url; do     curl -s "$url" | grep -o '<img[^>]*src="[^"]*"' | grep -o 'src="[^"]*"' | grep -o '"[^"]*"' | tr -d '"'; done
