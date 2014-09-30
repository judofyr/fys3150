set -e
pushd build

pdflatex ../article.tex

grep -q "Rerun to" article.log && pdflatex ../article.tex

if grep -q "Rerun to" article.log ; then
  echo What, again?!
  exit 1
fi

