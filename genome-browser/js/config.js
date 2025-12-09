/**
 * Lake Trout Genome Browser Configuration
 * IGV.js track and genome configuration
 */

// Base URL for data files - update this to your GitHub Pages URL
const DATA_BASE_URL = window.location.hostname === 'localhost' || window.location.hostname === '127.0.0.1'
    ? './data'
    : 'https://robertslab.github.io/project-lake-trout/genome-browser/data';

// External genome hosting URL (Gannet server)
const GENOME_BASE_URL = 'https://gannet.fish.washington.edu/v1_web/owlshell/bu-github/deep-dive-genome-browser/docs/jbrowse/data';

/**
 * Genome Configuration
 * Reference: GCF_016432855.1 SaNama_1.0 (Lake Trout / Salvelinus namaycush)
 * Note: Using full FASTA file for sequence display.
 */
const GENOME_CONFIG = {
    id: "SaNama_1.0",
    name: "Lake Trout (Salvelinus namaycush) - SaNama_1.0",
    // Full genome FASTA for sequence display
    fastaURL: `${GENOME_BASE_URL}/GCF_016432855.1_SaNama_1.0_genomic.fa`,
    indexURL: `${GENOME_BASE_URL}/GCF_016432855.1_SaNama_1.0_genomic.fa.fai`,
    // Explicitly set cytobandURL to null to prevent IGV.js errors
    cytobandURL: null,
    // Enable whole genome view with FASTA
    wholeGenomeView: true,
    
    // Chromosome aliases for user-friendly names
    chromosomeOrder: [
        "NC_052307.1", "NC_052308.1", "NC_052309.1", "NC_052310.1", "NC_052311.1",
        "NC_052312.1", "NC_052313.1", "NC_052314.1", "NC_052315.1", "NC_052316.1",
        "NC_052317.1", "NC_052318.1", "NC_052319.1", "NC_052320.1", "NC_052321.1",
        "NC_052322.1", "NC_052323.1", "NC_052324.1", "NC_052325.1", "NC_052326.1",
        "NC_052327.1", "NC_052328.1", "NC_052329.1", "NC_052330.1", "NC_052331.1",
        "NC_052332.1", "NC_052333.1", "NC_052334.1", "NC_052335.1", "NC_052336.1",
        "NC_052337.1", "NC_052338.1", "NC_052339.1", "NC_052340.1", "NC_052341.1",
        "NC_052342.1", "NC_052343.1", "NC_052344.1", "NC_052345.1", "NC_052346.1",
        "NC_052347.1", "NC_052348.1"
    ]
};

/**
 * Track Configurations
 * Organized by category for clear presentation
 */
const TRACK_CONFIGS = {
    // Gene Annotations
    genes: {
        name: "Genes",
        type: "annotation",
        format: "bed",
        url: `${DATA_BASE_URL}/annotations/genes.bed`,
        displayMode: "EXPANDED",
        color: "#10B981",
        height: 100,
        order: 1
    },
    
    // PAV Tracks - Lean-specific
    leanInsertions: {
        name: "Lean-Specific Insertions",
        type: "annotation",
        format: "bed",
        url: `${DATA_BASE_URL}/pav/lean_specific.insertions.bed`,
        displayMode: "SQUISHED",
        color: "#3B82F6",
        height: 50,
        order: 10,
        visibilityWindow: 10000000
    },
    leanDeletions: {
        name: "Lean-Specific Deletions", 
        type: "annotation",
        format: "bed",
        url: `${DATA_BASE_URL}/pav/lean_specific.deletions.bed`,
        displayMode: "SQUISHED",
        color: "#1D4ED8",
        height: 50,
        order: 11,
        visibilityWindow: 10000000
    },
    
    // PAV Tracks - Siscowet-specific
    siscowetInsertions: {
        name: "Siscowet-Specific Insertions",
        type: "annotation",
        format: "bed",
        url: `${DATA_BASE_URL}/pav/siscowet_specific.insertions.bed`,
        displayMode: "SQUISHED",
        color: "#EF4444",
        height: 50,
        order: 20,
        visibilityWindow: 10000000
    },
    siscowetDeletions: {
        name: "Siscowet-Specific Deletions",
        type: "annotation",
        format: "bed",
        url: `${DATA_BASE_URL}/pav/siscowet_specific.deletions.bed`,
        displayMode: "SQUISHED",
        color: "#B91C1C",
        height: 50,
        order: 21,
        visibilityWindow: 10000000
    },
    
    // Stringent PAV Tracks - Present in ALL 4 samples, absent from other ecotype, >100bp
    stringentSiscowetDeletions: {
        name: "Stringent Siscowet-Specific Deletions (>100bp, all 4 samples)",
        type: "annotation",
        format: "bed",
        url: `${DATA_BASE_URL}/pav/stringent.siscowet_specific.deletions.browser.bed`,
        displayMode: "SQUISHED",
        color: "#7F1D1D",
        height: 50,
        order: 22,
        visibilityWindow: 10000000
    },
    
    // PAV Tracks - Shared
    sharedInsertions: {
        name: "Shared Insertions",
        type: "annotation",
        format: "bed",
        url: `${DATA_BASE_URL}/pav/shared.insertions.bed`,
        displayMode: "SQUISHED",
        color: "#6B7280",
        height: 40,
        order: 30,
        visibilityWindow: 10000000
    },
    sharedDeletions: {
        name: "Shared Deletions",
        type: "annotation",
        format: "bed",
        url: `${DATA_BASE_URL}/pav/shared.deletions.bed`,
        displayMode: "SQUISHED",
        color: "#374151",
        height: 40,
        order: 31,
        visibilityWindow: 10000000
    },
    
    // Methylation Tracks
    methylationDiff: {
        name: "Methylation Difference (Siscowet - Lean)",
        type: "wig",
        format: "bedgraph",
        url: `${DATA_BASE_URL}/methylation/methylation_diff.bedGraph`,
        color: "#8B5CF6",
        altColor: "#F59E0B",
        height: 60,
        min: -50,
        max: 50,
        autoscale: false,
        order: 40,
        visibilityWindow: 5000000
    },
    
    // Individual sample methylation (optional - can be toggled)
    leanMethylation: {
        name: "Lean Mean Methylation",
        type: "wig",
        format: "bedgraph",
        url: `${DATA_BASE_URL}/methylation/lean_mean.bedGraph`,
        color: "#3B82F6",
        height: 50,
        min: 0,
        max: 100,
        order: 50,
        visibilityWindow: 2000000
    },
    siscowetMethylation: {
        name: "Siscowet Mean Methylation",
        type: "wig",
        format: "bedgraph",
        url: `${DATA_BASE_URL}/methylation/siscowet_mean.bedGraph`,
        color: "#EF4444",
        height: 50,
        min: 0,
        max: 100,
        order: 51,
        visibilityWindow: 2000000
    },
    
    // DMRs (Differentially Methylated Regions)
    dmrsHyper: {
        name: "DMRs - Hypermethylated in Siscowet",
        type: "annotation",
        format: "bed",
        url: `${DATA_BASE_URL}/methylation/dmrs_hyper_siscowet.bed`,
        displayMode: "SQUISHED",
        color: "#DC2626",
        height: 40,
        order: 60
    },
    dmrsHypo: {
        name: "DMRs - Hypomethylated in Siscowet",
        type: "annotation", 
        format: "bed",
        url: `${DATA_BASE_URL}/methylation/dmrs_hypo_siscowet.bed`,
        displayMode: "SQUISHED",
        color: "#2563EB",
        height: 40,
        order: 61
    },
    
    // Significant DMCs (Differentially Methylated CpGs)
    significantDMCs: {
        name: "Significant DMCs (p < 0.05)",
        type: "annotation",
        format: "bed",
        url: `${DATA_BASE_URL}/methylation/significant_dmcs.bed`,
        displayMode: "SQUISHED",
        color: "#A855F7",
        height: 40,
        order: 62,
        visibilityWindow: 5000000
    }
};

/**
 * Default track selection for initial load
 * Can be customized based on user preference
 */
const DEFAULT_TRACKS = [
    'genes',
    'leanInsertions',
    'leanDeletions', 
    'siscowetInsertions',
    'siscowetDeletions',
    'stringentSiscowetDeletions',
    'methylationDiff',
    'dmrsHyper',
    'dmrsHypo'
];

/**
 * Initial browser position
 */
const INITIAL_LOCUS = "NC_052307.1:1-500000";

/**
 * Browser options
 */
const BROWSER_OPTIONS = {
    showNavigation: true,
    showRuler: true,
    showCenterGuide: true,
    showCursorTrackingGuide: true,
    showSequence: false,  // Disabled - no FASTA file available
    flanking: 1000
};

// Export for use in app.js
window.GenomeBrowserConfig = {
    DATA_BASE_URL,
    GENOME_CONFIG,
    TRACK_CONFIGS,
    DEFAULT_TRACKS,
    INITIAL_LOCUS,
    BROWSER_OPTIONS
};

