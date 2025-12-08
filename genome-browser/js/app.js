/**
 * Lake Trout Genome Browser Application
 * Main application logic for IGV.js browser
 */

// Global browser instance
let browser = null;

/**
 * Initialize the genome browser
 */
async function initBrowser() {
    const config = window.GenomeBrowserConfig;
    
    // Show loading state
    const container = document.getElementById('igv-container');
    container.innerHTML = `
        <div class="loading-overlay">
            <div style="text-align: center;">
                <div class="loading-spinner"></div>
                <p style="margin-top: 1rem; color: var(--color-text-dim);">Loading genome browser...</p>
            </div>
        </div>
    `;
    
    // Check if config loaded properly
    if (!config || !config.GENOME_CONFIG) {
        container.innerHTML = `
            <div style="padding: 2rem; text-align: center; color: var(--color-text-dim);">
                <h3 style="color: var(--color-siscowet); margin-bottom: 1rem;">Error Loading Browser</h3>
                <p>Configuration failed to load.</p>
            </div>
        `;
        return;
    }
    
    try {
        // Build track list from default selection
        const tracks = config.DEFAULT_TRACKS
            .map(trackId => config.TRACK_CONFIGS[trackId])
            .filter(track => track !== undefined);
        
        // Determine if running locally or on GitHub Pages
        const isLocal = window.location.hostname === 'localhost' || window.location.hostname === '127.0.0.1';
        
        let referenceConfig;
        let showSequence;
        
        if (isLocal) {
            // Local development - use local FASTA files
            referenceConfig = {
                id: config.GENOME_CONFIG.id,
                name: config.GENOME_CONFIG.name,
                fastaURL: './data/genome/GCF_016432855.1_SaNama_1.0_genomic.fa',
                indexURL: './data/genome/GCF_016432855.1_SaNama_1.0_genomic.fa.fai'
            };
            showSequence = true;
        } else {
            // GitHub Pages - define chromosomes inline to avoid external file loading issues
            // This avoids the cytobandURL error with chromSizesURL
            referenceConfig = {
                id: config.GENOME_CONFIG.id,
                name: config.GENOME_CONFIG.name,
                chromosomes: [
                    {name: "NC_052307.1", size: 84126519},
                    {name: "NC_052308.1", size: 80899063},
                    {name: "NC_052309.1", size: 98200354},
                    {name: "NC_052310.1", size: 84641001},
                    {name: "NC_052311.1", size: 74663310},
                    {name: "NC_052312.1", size: 53240245},
                    {name: "NC_052313.1", size: 62955435},
                    {name: "NC_052314.1", size: 89946781},
                    {name: "NC_052315.1", size: 77627323},
                    {name: "NC_052316.1", size: 55315881},
                    {name: "NC_052317.1", size: 72888909},
                    {name: "NC_052318.1", size: 76979435},
                    {name: "NC_052319.1", size: 78589171},
                    {name: "NC_052320.1", size: 64649498},
                    {name: "NC_052321.1", size: 65048267},
                    {name: "NC_052322.1", size: 55422761},
                    {name: "NC_052323.1", size: 52735498},
                    {name: "NC_052324.1", size: 48040928},
                    {name: "NC_052325.1", size: 56779286},
                    {name: "NC_052326.1", size: 42509555},
                    {name: "NC_052327.1", size: 49837422},
                    {name: "NC_052328.1", size: 50166449},
                    {name: "NC_052329.1", size: 48355507},
                    {name: "NC_052330.1", size: 47124051},
                    {name: "NC_052331.1", size: 42984373},
                    {name: "NC_052332.1", size: 52270378},
                    {name: "NC_052333.1", size: 40490003},
                    {name: "NC_052334.1", size: 42540729},
                    {name: "NC_052335.1", size: 50269628},
                    {name: "NC_052336.1", size: 44632854},
                    {name: "NC_052337.1", size: 41213638},
                    {name: "NC_052338.1", size: 41426868},
                    {name: "NC_052339.1", size: 40758918},
                    {name: "NC_052340.1", size: 31687327},
                    {name: "NC_052341.1", size: 35142098},
                    {name: "NC_052342.1", size: 31730470},
                    {name: "NC_052343.1", size: 36544277},
                    {name: "NC_052344.1", size: 42306528},
                    {name: "NC_052345.1", size: 31458917},
                    {name: "NC_052346.1", size: 37618679},
                    {name: "NC_052347.1", size: 36174209},
                    {name: "NC_052348.1", size: 28563236}
                ]
            };
            showSequence = false;
        }
        
        // IGV.js configuration
        const igvConfig = {
            showNavigation: true,
            showRuler: true,
            showCenterGuide: true,
            showCursorTrackingGuide: true,
            showSequence: showSequence,
            locus: config.INITIAL_LOCUS,
            reference: referenceConfig,
            tracks: tracks
        };
        
        // Debug logging
        console.log('IGV Config:', JSON.stringify(igvConfig, null, 2));
        
        // Clear container and create browser
        container.innerHTML = '';
        browser = await igv.createBrowser(container, igvConfig);
        
        console.log('Lake Trout Genome Browser initialized successfully');
        
        // Set up event handlers
        setupEventHandlers();
        
    } catch (error) {
        console.error('Error initializing browser:', error);
        container.innerHTML = `
            <div style="padding: 2rem; text-align: center; color: var(--color-text-dim);">
                <h3 style="color: var(--color-siscowet); margin-bottom: 1rem;">Error Loading Browser</h3>
                <p>${error.message}</p>
                <p style="margin-top: 1rem; font-size: 0.875rem;">
                    This may occur if data files are not yet prepared.<br>
                    Run <code style="background: var(--color-bg-elevated); padding: 2px 6px; border-radius: 4px;">python prepare_data.py</code> 
                    to generate track files.
                </p>
            </div>
        `;
    }
}

/**
 * Set up event handlers for UI interactions
 */
function setupEventHandlers() {
    // Region search input
    const searchInput = document.getElementById('region-search');
    if (searchInput) {
        searchInput.addEventListener('keypress', (e) => {
            if (e.key === 'Enter') {
                const region = searchInput.value.trim();
                if (region) {
                    navigateToRegion(region);
                }
            }
        });
    }
    
    // Browser track change events
    if (browser) {
        browser.on('locuschange', (referenceFrame) => {
            // Update URL hash for bookmarking
            if (referenceFrame && referenceFrame.length > 0) {
                const locus = referenceFrame[0].getLocusString();
                window.history.replaceState(null, '', `#${locus}`);
            }
        });
    }
    
    // Handle URL hash on load
    handleUrlHash();
}

/**
 * Navigate to a specific genomic region
 * @param {string} region - Region string (e.g., "NC_052307.1:100000-200000")
 */
function navigateToRegion(region) {
    if (browser) {
        try {
            browser.search(region);
        } catch (error) {
            console.error('Error navigating to region:', error);
            alert(`Could not navigate to region: ${region}`);
        }
    }
}

/**
 * Handle URL hash for deep linking
 */
function handleUrlHash() {
    const hash = window.location.hash.substring(1);
    if (hash && browser) {
        setTimeout(() => {
            navigateToRegion(hash);
        }, 500);
    }
}

/**
 * Add a track to the browser
 * @param {string} trackId - Track identifier from config
 */
async function addTrack(trackId) {
    const config = window.GenomeBrowserConfig;
    const trackConfig = config.TRACK_CONFIGS[trackId];
    
    if (!trackConfig) {
        console.error(`Track not found: ${trackId}`);
        return;
    }
    
    if (browser) {
        try {
            await browser.loadTrack(trackConfig);
            console.log(`Track added: ${trackConfig.name}`);
        } catch (error) {
            console.error(`Error adding track ${trackId}:`, error);
        }
    }
}

/**
 * Remove a track from the browser
 * @param {string} trackName - Name of the track to remove
 */
function removeTrack(trackName) {
    if (browser) {
        const tracks = browser.trackViews;
        for (const trackView of tracks) {
            if (trackView.track.name === trackName) {
                browser.removeTrack(trackView.track);
                console.log(`Track removed: ${trackName}`);
                return;
            }
        }
    }
}

/**
 * Get all available tracks
 */
function getAvailableTracks() {
    return window.GenomeBrowserConfig.TRACK_CONFIGS;
}

/**
 * Export current view as SVG
 */
async function exportSVG() {
    if (browser) {
        try {
            const svg = await browser.toSVG();
            const blob = new Blob([svg], { type: 'image/svg+xml' });
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = 'lake_trout_browser.svg';
            a.click();
            URL.revokeObjectURL(url);
        } catch (error) {
            console.error('Error exporting SVG:', error);
        }
    }
}

/**
 * Utility: Format genomic position for display
 * @param {number} position - Genomic position
 */
function formatPosition(position) {
    return position.toLocaleString();
}

/**
 * Utility: Parse region string
 * @param {string} region - Region string (e.g., "chr1:100-200")
 */
function parseRegion(region) {
    const match = region.match(/^(.+):(\d+)-(\d+)$/);
    if (match) {
        return {
            chromosome: match[1],
            start: parseInt(match[2]),
            end: parseInt(match[3])
        };
    }
    return null;
}

// Initialize browser when DOM is ready
document.addEventListener('DOMContentLoaded', initBrowser);

// Make functions globally available
window.navigateToRegion = navigateToRegion;
window.addTrack = addTrack;
window.removeTrack = removeTrack;
window.getAvailableTracks = getAvailableTracks;
window.exportSVG = exportSVG;

