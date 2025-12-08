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
        
        // Gannet URL for genome files (repo syncs to Gannet)
        const gannetGenomeURL = 'https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/project-lake-trout/genome-browser/data/genome';
        
        let referenceConfig;
        let showSequence;
        
        if (isLocal) {
            // Local development - use local files
            referenceConfig = {
                id: config.GENOME_CONFIG.id,
                name: config.GENOME_CONFIG.name,
                fastaURL: './data/genome/GCF_016432855.1_SaNama_1.0_genomic.fa',
                indexURL: './data/genome/GCF_016432855.1_SaNama_1.0_genomic.fa.fai'
            };
            showSequence = true;
        } else {
            // GitHub Pages - use Gannet for genome files (repo syncs there)
            referenceConfig = {
                id: config.GENOME_CONFIG.id,
                name: config.GENOME_CONFIG.name,
                fastaURL: gannetGenomeURL + '/GCF_016432855.1_SaNama_1.0_genomic.fa',
                indexURL: gannetGenomeURL + '/GCF_016432855.1_SaNama_1.0_genomic.fa.fai'
            };
            showSequence = true;
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

