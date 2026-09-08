#!/usr/bin/env bash
set -euo pipefail

: "${CTAN_URL:=https://mirrors.rit.edu/CTAN}"
: "${TL_PACKAGES:=accsupp amscls amsmath anyfontsize blkarray booktabs braket caption chemformula ctex endnotes enumitem fancyhdr float graphics hyperref latexmk luatex85 l3kernel l3packages l3experimental mathtools metafont multirow needspace newfloat pgfplots scalerel siunitx stackengine threeparttable tools ulem xcolor xecjk xfrac}"
: "${TL_FONT_PCK:=boondox fandol libertinus-fonts mathalpha psnfss collection-fontsrecommended}"

export PATH=/tmp/texlive/bin/x86_64-linux:$PATH

TLREPO="$CTAN_URL/systems/texlive/tlnet"
PROFILE="${GITHUB_WORKSPACE:-$PWD}/.github/workflows/texlive.profile"

wget "$TLREPO/install-tl-unx.tar.gz"
tar -xzf install-tl-unx.tar.gz
cd install-tl-20*
./install-tl --profile "$PROFILE"
tlmgr option repository "$TLREPO"
tlmgr install $TL_PACKAGES
tlmgr install $TL_FONT_PCK
tlmgr update --self --all --no-auto-install --repository="$TLREPO"
