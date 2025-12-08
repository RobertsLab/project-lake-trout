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
    
    try {
        // Build track list from default selection
        const tracks = config.DEFAULT_TRACKS
            .map(trackId => config.TRACK_CONFIGS[trackId])
            .filter(track => track !== undefined);
        
        // Build reference config - use FASTA if available, otherwise chromSizes
        const referenceConfig = {
            id: config.GENOME_CONFIG.id,
            name: config.GENOME_CONFIG.name,
            chromosomeOrder: config.GENOME_CONFIG.chromosomeOrder
        };
        
        // Add FASTA or chromSizes based on what's configured
        if (config.GENOME_CONFIG.fastaURL) {
            referenceConfig.fastaURL = config.GENOME_CONFIG.fastaURL;
            referenceConfig.indexURL = config.GENOME_CONFIG.indexURL;
        } else if (config.GENOME_CONFIG.chromSizesURL) {
            referenceConfig.chromSizesURL = config.GENOME_CONFIG.chromSizesURL;
        }
        
        // IGV.js configuration
        const igvConfig = {
            ...config.BROWSER_OPTIONS,
            locus: config.INITIAL_LOCUS,
            reference: referenceConfig,
            tracks: tracks
        };
        
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

