

# Created with

  docker run --rm --volume=$(pwd):/srv/jekyll -i -t  -p 127.0.0.1:4000:4000 jekyll/jekyll:pages jekyll new --force .

# Local hosting

```
docker run --rm --volume=$(pwd):/srv/jekyll -i -t  -p 127.0.0.1:4000:4000 jekyll/jekyll:pages
```

The website will be hosted at `http://localhost:4000/<repo name>/`.
