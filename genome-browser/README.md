# Lake Trout Genome Browser

Interactive genome browser for visualizing **Presence-Absence Variation (PAV)** and **DNA methylation** differences between Lake Trout (*Salvelinus namaycush*) ecotypes: **Lean** and **Siscowet**.

## Features

- **PAV Visualization**: View insertions and deletions specific to each ecotype
- **Methylation Tracks**: Visualize differential methylation between ecotypes
- **Gene Annotations**: Reference gene annotations from SaNama_1.0 assembly
- **Interactive Navigation**: Jump to any genomic region
- **Mobile Responsive**: Works on desktop and mobile devices

## Live Demo

Visit: **https://sr320.github.io/project-lake-trout/genome-browser/**

## Quick Start

### 1. Prepare Data

First, run the prerequisite analyses (if not already done):

```bash
# Generate differential PAV data
cd /path/to/project-lake-trout
python code/15-diff-pav.py

# Generate differential methylation data (optional)
python code/14-diff-meth.py
```

Then prepare the browser data files:

```bash
cd genome-browser
python prepare_data.py
```

### 2. Local Development

Start a local server:

```bash
cd genome-browser
python -m http.server 8000
```

Open http://localhost:8000 in your browser.

### 3. Deploy to GitHub Pages

1. Commit and push changes:
```bash
git add genome-browser/
git commit -m "Add genome browser"
git push origin main
```

2. Enable GitHub Pages:
   - Go to repository Settings → Pages
   - Source: Deploy from branch `main`
   - Folder: `/ (root)` or `/docs` if moved there

3. Access at: `https://[username].github.io/project-lake-trout/genome-browser/`

## Directory Structure

```
genome-browser/
├── index.html              # Main HTML page
├── css/
│   └── style.css           # Stylesheet
├── js/
│   ├── config.js           # Track configurations
│   └── app.js              # Application logic
├── data/                   # Prepared data files (generated)
│   ├── genome/             # Reference genome (symlinks)
│   ├── annotations/        # Gene annotations
│   ├── pav/                # PAV BED files
│   └── methylation/        # Methylation bedGraph files
├── prepare_data.py         # Data preparation script
└── README.md               # This file
```

## Data Files

### PAV Tracks

| Track | Description | Count | Color |
|-------|-------------|-------|-------|
| Lean-Specific Insertions | Insertions found only in Lean | 770,891 | Blue |
| Lean-Specific Deletions | Deletions found only in Lean | 225,337 | Dark Blue |
| Siscowet-Specific Insertions | Insertions found only in Siscowet | 1,086,799 | Red |
| Siscowet-Specific Deletions | Deletions found only in Siscowet | 245,906 | Dark Red |
| Shared Insertions | Present in both ecotypes | 672,773 | Gray |
| Shared Deletions | Present in both ecotypes | 205,599 | Dark Gray |

### Methylation Tracks

| Track | Description | Count |
|-------|-------------|-------|
| Methylation Difference | Siscowet - Lean methylation % (positive = higher in Siscowet) | 540,040 CpGs tested |
| DMRs Hypermethylated | Regions with higher methylation in Siscowet | 20 |
| DMRs Hypomethylated | Regions with lower methylation in Siscowet | 282 |
| Significant DMCs | Individual CpG sites with significant differential methylation (p < 0.05) | 4,440 |

### Analysis Summary

**Differential PAV Analysis:**
- Lean-specific: 996,228 variants (770,891 insertions + 225,337 deletions)
- Siscowet-specific: 1,332,705 variants (1,086,799 insertions + 245,906 deletions)
- Shared between ecotypes: 878,372 variants (672,773 insertions + 205,599 deletions)

**Differential Methylation Analysis:**
- CpG sites tested: 540,040
- Significant DMCs (p < 0.05): 4,440
  - Hypermethylated in Siscowet: 445
  - Hypomethylated in Siscowet: 3,995
- Total DMRs: 302
  - Hypermethylated in Siscowet: 20
  - Hypomethylated in Siscowet: 282
- Mean methylation difference: 23.2%

## Sample Information

### Lean Ecotype
- bc2041, bc2069, bc2070, bc2068

### Siscowet Ecotype
- bc2071, bc2073, bc2072, bc2096

## Reference Genome

- **Assembly**: GCF_016432855.1
- **Name**: SaNama_1.0
- **Species**: *Salvelinus namaycush* (Lake Trout)

## Hosting Large Files

For GitHub Pages deployment, large genome files (FASTA) may need external hosting:

### Option 1: Use NCBI (recommended)
Configure IGV.js to load from NCBI servers directly.

### Option 2: External Server
Host on a CORS-enabled web server and update `DATA_BASE_URL` in `js/config.js`.

### Option 3: Gannet Server
If you have access to the Gannet server:
```javascript
const DATA_BASE_URL = 'https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/project-lake-trout/genome-browser/data';
```

## Technology Stack

- **[IGV.js](https://github.com/igvteam/igv.js)** - Interactive genome visualization
- **Pure CSS** - No framework dependencies
- **Vanilla JavaScript** - No build step required

## Customization

### Adding New Tracks

Edit `js/config.js` and add to `TRACK_CONFIGS`:

```javascript
myNewTrack: {
    name: "My New Track",
    type: "annotation",
    format: "bed",
    url: `${DATA_BASE_URL}/mydata/track.bed`,
    displayMode: "EXPANDED",
    color: "#FF0000",
    height: 50,
    order: 100
}
```

Then add to `DEFAULT_TRACKS` array to show by default.

### Changing Colors

Edit the CSS variables in `css/style.css`:

```css
:root {
    --color-lean: #3B82F6;      /* Blue */
    --color-siscowet: #EF4444;  /* Red */
    --color-accent: #06B6D4;    /* Cyan */
}
```

## Troubleshooting

### Browser shows "Error Loading Browser"
- Run `python prepare_data.py` to generate data files
- Check browser console for specific error messages
- Verify data files exist in `data/` directory

### Tracks not loading
- Check CORS settings if hosting externally
- Verify file URLs in `js/config.js`
- Ensure BED files have correct format

### Slow performance
- Large BED files may need bigBed conversion
- Increase `visibilityWindow` to limit features shown
- Use track `displayMode: "COLLAPSED"` for dense tracks

## License

This project is part of the [project-lake-trout](https://github.com/sr320/project-lake-trout) repository.

## Citation

If you use this browser in your research, please cite:
- The project-lake-trout repository
- IGV.js: Robinson et al., 2023

## Contact

For questions or issues, please open a GitHub issue in the main repository.

