#!/usr/bin/env python3
"""
ParaChrSNP Web Launcher

Author: Junpeng Ma 1527552938@qq.com
"""

import argparse
import gzip
import html
import json
import mimetypes
import os
import re
import shlex
import signal
import subprocess
import sys
import threading
import time
import uuid
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

try:
    import yaml
except ImportError:
    yaml = None


AUTHOR = "Junpeng Ma 1527552938@qq.com"
SCRIPT_NAME = "parachrsnp_web.py"
SCRIPT_FUNCTION = "Provide a browser-based launcher for ParaChrSNP Snakemake workflow."
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = PROJECT_ROOT / "web_runs"
DEFAULT_SNAKEMAKE = str(Path(sys.executable).resolve().parent / "snakemake")
DEFAULT_ALLOWED_ROOTS = [
    PROJECT_ROOT,
    Path("/home/majunpeng/sda2"),
    Path("/data"),
    Path("/public"),
]

JOBS = {}
JOBS_LOCK = threading.Lock()

STAGE_RULES = {
    "precheck": ["precheck"],
    "qc": ["qc"],
    "clean_reads": ["clean_reads"],
    "bwa_mem": ["bwa_index", "bwa_mem", "faidx_index", "picard_index"],
    "bam_rmdup": ["bam_rmdup", "index_rmdup"],
    "haplo": ["haplo", "get_chr_list", "merge_sample_gvcf", "combine_gvcf"],
    "joint_calling": ["genomicsdb_sample_map", "genomicsdb_import", "genotype_chrom", "joint_calling"],
    "snp_filter": ["snp_select", "snp_filter", "indel_select", "indel_filter"],
    "vcf_missing": ["vcf_missing"],
    "vcf_to": ["vcf_to_plink_binary", "vcf_to_plink_text", "vcf_to_hapmap"],
    "imputation": ["imputation"],
    "cnv": ["cnvnator_call", "cnv_filter", "cnv_summary"],
    "vcf2pca": ["vcf2pca"],
    "vcf2dis": ["vcf2dis"],
    "snpeff": ["snpeff_prepare_database", "snpeff_build_database", "snpeff_annotate_snp", "snpeff_annotate_indel"],
    "report": ["report"],
}


INDEX_HTML = r"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ParaChrSNP Web Launcher</title>
  <style>
    :root {
      --bg: #edf1f4;
      --panel: #ffffff;
      --panel2: #f8fafb;
      --text: #1f2933;
      --muted: #687782;
      --line: #d8e0e5;
      --accent: #0e756b;
      --accent-dark: #095c55;
      --accent-soft: #e3f1ef;
      --warn: #b65a31;
      --danger: #b83f3f;
      --ink: #162129;
      --shadow: 0 10px 30px rgba(25, 43, 58, .08);
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      background: var(--bg);
      color: var(--text);
      font-family: Arial, "Noto Sans SC", "Microsoft YaHei", sans-serif;
      font-size: 14px;
      line-height: 1.5;
    }
    .app-shell {
      min-height: 100vh;
      display: grid;
      grid-template-rows: auto 1fr;
    }
    header {
      background: rgba(255,255,255,.94);
      border-bottom: 1px solid var(--line);
      padding: 14px 26px;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 16px;
      position: sticky;
      top: 0;
      z-index: 10;
      backdrop-filter: blur(10px);
    }
    .brand {
      display: flex;
      align-items: center;
      gap: 12px;
      min-width: 220px;
    }
    .brand img {
      width: 42px;
      height: 42px;
      object-fit: contain;
    }
    .brand h1 {
      margin: 0;
      font-size: 22px;
      font-weight: 700;
      letter-spacing: 0;
    }
    .subtitle {
      color: var(--muted);
      font-size: 13px;
      margin-top: 1px;
    }
    .status-bar {
      display: flex;
      gap: 10px;
      align-items: center;
      flex-wrap: wrap;
      justify-content: flex-end;
    }
    .badge {
      border: 1px solid var(--line);
      background: var(--panel2);
      color: var(--ink);
      padding: 5px 10px;
      border-radius: 6px;
      font-size: 13px;
    }
    main {
      display: grid;
      grid-template-columns: minmax(390px, 540px) minmax(480px, 1fr);
      gap: 18px;
      padding: 18px 26px 26px;
    }
    section {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 18px;
      min-width: 0;
      box-shadow: var(--shadow);
    }
    h2 {
      margin: 0 0 12px;
      font-size: 16px;
      font-weight: 700;
      letter-spacing: 0;
    }
    h3 {
      margin: 16px 0 8px;
      font-size: 14px;
      font-weight: 700;
    }
    .section-head {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      margin-bottom: 12px;
    }
    .section-head h2 { margin: 0; }
    label {
      display: block;
      font-weight: 600;
      margin: 10px 0 5px;
    }
    input, textarea, select {
      width: 100%;
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 8px 10px;
      background: #fff;
      color: var(--text);
      font: inherit;
      min-height: 36px;
    }
    input:focus, textarea:focus, select:focus {
      outline: 2px solid rgba(14,117,107,.16);
      border-color: var(--accent);
    }
    textarea {
      min-height: 96px;
      resize: vertical;
      font-family: "DejaVu Sans Mono", Consolas, monospace;
      font-size: 13px;
    }
    .row {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 10px;
    }
    .path-row {
      display: grid;
      grid-template-columns: 1fr auto;
      gap: 8px;
      align-items: end;
    }
    .actions {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      margin-top: 14px;
    }
    button {
      border: 1px solid var(--accent);
      background: var(--accent);
      color: #fff;
      border-radius: 6px;
      min-height: 36px;
      padding: 7px 12px;
      cursor: pointer;
      font-weight: 600;
    }
    button:hover { background: var(--accent-dark); }
    button.secondary {
      background: #fff;
      color: var(--accent);
    }
    button.secondary:hover { background: var(--accent-soft); }
    button.ghost {
      background: var(--panel2);
      border-color: var(--line);
      color: var(--ink);
    }
    button.ghost:hover { background: #eef4f3; }
    button.danger {
      border-color: var(--danger);
      background: var(--danger);
    }
    button:disabled {
      cursor: not-allowed;
      opacity: .55;
    }
    .browse-btn {
      min-width: 74px;
      padding-inline: 10px;
    }
    .checks {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 8px;
      margin-top: 8px;
    }
    .check {
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 8px 10px;
      display: flex;
      align-items: center;
      gap: 8px;
      min-height: 38px;
      background: var(--panel2);
    }
    .check input {
      width: auto;
      min-height: 0;
    }
    .advanced-panel {
      margin-top: 14px;
      background: #fbfdfd;
    }
    .param-grid {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 10px;
      margin-top: 8px;
    }
    .param-grid label {
      margin-top: 0;
    }
    .param-note {
      color: var(--muted);
      font-size: 12px;
      margin: 6px 0 0;
    }
    .dashboard {
      display: grid;
      grid-template-columns: repeat(4, minmax(110px, 1fr));
      gap: 10px;
      margin-bottom: 14px;
    }
    .metric {
      border: 1px solid var(--line);
      background: var(--panel2);
      border-radius: 8px;
      padding: 11px;
      min-height: 72px;
    }
    .metric .value {
      font-size: 22px;
      font-weight: 700;
      color: var(--ink);
      line-height: 1.1;
    }
    .metric .label {
      color: var(--muted);
      font-size: 12px;
      margin-top: 4px;
    }
    .progress-wrap {
      border: 1px solid var(--line);
      background: #fff;
      border-radius: 8px;
      padding: 14px;
      margin-bottom: 14px;
    }
    .report-actions {
      display: none;
      border: 1px solid var(--line);
      background: #f7fbfa;
      border-radius: 8px;
      padding: 12px;
      margin-bottom: 14px;
      align-items: center;
      justify-content: space-between;
      gap: 10px;
    }
    .report-actions.visible { display: flex; }
    .report-actions .muted { margin: 0; }
    .progress-header {
      display: flex;
      justify-content: space-between;
      gap: 12px;
      align-items: center;
      margin-bottom: 8px;
    }
    .progress-bar {
      height: 14px;
      border-radius: 999px;
      background: #e6ecef;
      overflow: hidden;
      border: 1px solid #d7e1e5;
    }
    .progress-fill {
      height: 100%;
      width: 0%;
      background: linear-gradient(90deg, #0e756b, #3a9c8a);
      transition: width .25s ease;
    }
    .stage-list {
      display: grid;
      grid-template-columns: repeat(3, 1fr);
      gap: 8px;
      margin-top: 10px;
    }
    .stage {
      border: 2px solid #2f7ecb;
      outline: 1px solid #2f7ecb;
      outline-offset: -4px;
      background: var(--panel2);
      border-radius: 8px;
      padding: 8px 10px;
      min-height: 48px;
      transition: border-color .2s ease, outline-color .2s ease, background .2s ease, box-shadow .2s ease;
    }
    .stage .name {
      font-weight: 700;
      font-size: 13px;
    }
    .stage .state {
      color: var(--muted);
      font-size: 12px;
      margin-top: 2px;
    }
    .stage.pending {
      border-color: #2f7ecb;
      outline-color: #2f7ecb;
      background: #f4f9ff;
    }
    .stage.completed {
      border-color: #25965f;
      outline-color: #25965f;
      background: #f2fbf6;
    }
    .stage.running {
      border-color: #ca3b3b;
      outline-color: #ca3b3b;
      background: #fff5f5;
      animation: runningPulse 1.1s ease-in-out infinite;
    }
    .stage.completed .state { color: #1f7a4d; }
    .stage.running .state { color: #b42f2f; font-weight: 700; }
    .stage.pending .state { color: #2f6fae; }
    @keyframes runningPulse {
      0%, 100% { box-shadow: 0 0 0 0 rgba(202,59,59,.18); }
      50% { box-shadow: 0 0 0 4px rgba(202,59,59,.16); }
    }
    table {
      width: 100%;
      border-collapse: collapse;
      font-size: 13px;
      margin-top: 8px;
    }
    th, td {
      border-bottom: 1px solid var(--line);
      padding: 7px 6px;
      text-align: left;
      vertical-align: top;
      word-break: break-all;
    }
    th {
      background: #f7faf9;
      font-weight: 700;
    }
    .table-wrap {
      max-height: 220px;
      overflow: auto;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: #fff;
    }
    .table-wrap table { margin-top: 0; }
    .log {
      background: #111820;
      color: #dbe7e3;
      border-radius: 8px;
      padding: 12px;
      min-height: 260px;
      max-height: 46vh;
      overflow: auto;
      white-space: pre-wrap;
      font-family: "DejaVu Sans Mono", Consolas, monospace;
      font-size: 12px;
    }
    .muted {
      color: var(--muted);
      font-size: 13px;
    }
    .split {
      display: grid;
      grid-template-rows: auto 1fr;
      gap: 12px;
    }
    details {
      border: 1px solid var(--line);
      border-radius: 8px;
      background: #fff;
      padding: 10px 12px;
    }
    summary {
      cursor: pointer;
      font-weight: 700;
      color: var(--ink);
    }
    .modal {
      position: fixed;
      inset: 0;
      background: rgba(18, 29, 38, .44);
      display: none;
      align-items: center;
      justify-content: center;
      z-index: 50;
      padding: 18px;
    }
    .modal.open { display: flex; }
    .browser {
      background: #fff;
      border-radius: 10px;
      border: 1px solid var(--line);
      box-shadow: 0 22px 70px rgba(0,0,0,.24);
      width: min(920px, 96vw);
      max-height: 86vh;
      display: grid;
      grid-template-rows: auto auto 1fr auto;
      overflow: hidden;
    }
    .browser-head {
      padding: 14px 16px;
      border-bottom: 1px solid var(--line);
      display: flex;
      justify-content: space-between;
      gap: 12px;
      align-items: center;
    }
    .browser-head h2 { margin: 0; }
    .browser-path {
      padding: 10px 16px;
      background: var(--panel2);
      border-bottom: 1px solid var(--line);
      display: grid;
      grid-template-columns: 1fr auto;
      gap: 8px;
      align-items: center;
    }
    .browser-list {
      overflow: auto;
      min-height: 280px;
      max-height: 52vh;
    }
    .file-row {
      display: grid;
      grid-template-columns: 92px 1fr 112px 160px;
      gap: 10px;
      padding: 9px 16px;
      border-bottom: 1px solid var(--line);
      align-items: center;
      cursor: pointer;
    }
    .file-row:hover { background: var(--accent-soft); }
    .file-row.selected {
      background: #d9eeeb;
      outline: 1px solid rgba(14,117,107,.36);
    }
    .file-kind {
      font-weight: 700;
      color: var(--accent-dark);
    }
    .file-name { word-break: break-all; }
    .browser-foot {
      padding: 12px 16px;
      border-top: 1px solid var(--line);
      display: flex;
      justify-content: space-between;
      gap: 12px;
      align-items: center;
    }
    @media (max-width: 980px) {
      main { grid-template-columns: 1fr; padding: 12px; }
      header { align-items: flex-start; flex-direction: column; }
      .status-bar { justify-content: flex-start; }
      .row, .checks, .dashboard, .stage-list { grid-template-columns: 1fr; }
      .file-row { grid-template-columns: 72px 1fr; }
      .file-size, .file-time { display: none; }
    }
  </style>
</head>
<body>
  <div class="app-shell">
    <header>
      <div class="brand">
        <img src="/icon.png" alt="ParaChrSNP">
        <div>
          <h1>ParaChrSNP</h1>
          <div class="subtitle">服务器端可视化工作流控制台</div>
        </div>
      </div>
      <div class="status-bar">
        <span class="badge" id="serverStatus">server: checking</span>
        <span class="badge" id="jobStatus">job: idle</span>
      </div>
    </header>

    <main>
      <section>
        <div class="section-head">
          <h2>输入与参数</h2>
          <span class="badge">server paths</span>
        </div>

        <label>FASTQ 目录</label>
        <div class="path-row">
          <input id="fastqDir" value="raw_fastq">
          <button class="secondary browse-btn" onclick="openBrowser('fastqDir','dir')">浏览</button>
        </div>
        <div class="actions">
          <button class="ghost" onclick="inferSamples()">识别样本</button>
        </div>

        <label>参考基因组 FASTA</label>
        <div class="path-row">
          <input id="reference" value="reference/Arabidopsis_thaliana.fasta">
          <button class="secondary browse-btn" onclick="openBrowser('reference','file','.fa,.fasta,.fna,.fa.gz,.fasta.gz')">浏览</button>
        </div>
        <div class="actions">
          <button class="ghost" onclick="inferChromosomes()">读取染色体</button>
        </div>

        <label>GFF/GTF 注释文件</label>
        <div class="path-row">
          <input id="annotationFile" value="reference/Arabidopsis_thaliana.gff">
          <button class="secondary browse-btn" onclick="openBrowser('annotationFile','file','.gff,.gff3,.gtf,.gff.gz,.gff3.gz,.gtf.gz')">浏览</button>
        </div>

        <div class="row">
          <div>
            <label>容器镜像</label>
            <div class="path-row">
              <input id="containerImage" value="ParaChrSNP_ape_ggtree.sif">
              <button class="secondary browse-btn" onclick="openBrowser('containerImage','file','.sif,.simg')">浏览</button>
            </div>
          </div>
          <div>
            <label>CPU 核心数</label>
            <input id="cores" type="number" min="1" value="12">
          </div>
        </div>

        <label>Snakemake 命令</label>
        <input id="snakemakeCmd" value="__DEFAULT_SNAKEMAKE__">

        <label>Singularity 绑定参数</label>
        <input id="singularityArgs" value="-B /home/majunpeng/sda2">

        <h3>样本</h3>
        <textarea id="samples" placeholder="每行一个样本：sample_name<TAB>fastq_prefix"></textarea>

        <h3>染色体</h3>
        <textarea id="chromosomes" placeholder="每行一个染色体 ID"></textarea>

        <h3>可选模块</h3>
        <div class="checks">
          <label class="check"><input id="runPca" type="checkbox" checked> PCA</label>
          <label class="check"><input id="runDis" type="checkbox" checked> 遗传距离/树</label>
          <label class="check"><input id="runSnpeff" type="checkbox"> SnpEff 注释</label>
          <label class="check"><input id="runIndelAnno" type="checkbox"> 注释 INDEL</label>
          <label class="check"><input id="runImputation" type="checkbox"> Beagle 填充</label>
          <label class="check"><input id="runCnv" type="checkbox"> CNVnator CNV</label>
        </div>

        <details class="advanced-panel" open>
          <summary>高级参数</summary>

          <h3>线程数</h3>
          <div class="param-grid">
            <div>
              <label>RastQC 线程数</label>
              <input id="qcThreads" type="number" min="1" value="8">
            </div>
            <div>
              <label>fastp 线程数</label>
              <input id="fastpThreads" type="number" min="1" value="8">
            </div>
            <div>
              <label>bwa-mem2 比对线程数</label>
              <input id="bwaThreads" type="number" min="1" value="8">
            </div>
            <div>
              <label>VCF2PCACluster 线程数</label>
              <input id="pcaThreads" type="number" min="1" value="4">
            </div>
            <div>
              <label>Beagle 线程数</label>
              <input id="imputationThreads" type="number" min="1" value="4">
            </div>
            <div>
              <label>GenomicsDBImport 读取线程数</label>
              <input id="genomicsdbThreads" type="number" min="1" value="4">
            </div>
            <div>
              <label>GenomicsDBImport batch size</label>
              <input id="genomicsdbBatchSize" type="number" min="1" value="50">
            </div>
            <div>
              <label>CNVnator bin size</label>
              <input id="cnvBinSize" type="number" min="1" value="100">
            </div>
            <div>
              <label>CNV 最小长度</label>
              <input id="cnvMinLen" type="number" min="1" value="1000">
            </div>
            <div>
              <label>CNV 最大 e-value</label>
              <input id="cnvMaxEvalue" value="0.01">
            </div>
            <div>
              <label>CNV 最大 q0</label>
              <input id="cnvMaxQ0" value="0.5">
            </div>
          </div>

          <h3>Java 内存</h3>
          <div class="param-grid">
            <div>
              <label>HaplotypeCaller 最小内存</label>
              <input id="haploXms" value="512m">
            </div>
            <div>
              <label>HaplotypeCaller 最大内存</label>
              <input id="haploXmx" value="4g">
            </div>
            <div>
              <label>GenotypeGVCFs 最小内存</label>
              <input id="genotypeXms" value="512m">
            </div>
            <div>
              <label>GenotypeGVCFs 最大内存</label>
              <input id="genotypeXmx" value="128g">
            </div>
            <div>
              <label>CombineGVCFs 最小内存</label>
              <input id="combineXms" value="512m">
            </div>
            <div>
              <label>CombineGVCFs 最大内存</label>
              <input id="combineXmx" value="128g">
            </div>
            <div>
              <label>GenomicsDBImport 最小内存</label>
              <input id="genomicsdbXms" value="1g">
            </div>
            <div>
              <label>GenomicsDBImport 最大内存</label>
              <input id="genomicsdbXmx" value="16g">
            </div>
            <div>
              <label>SnpEff 最大内存</label>
              <input id="snpeffXmx" value="32g">
            </div>
            <div>
              <label>Beagle 最大内存</label>
              <input id="beagleXmx" value="4g">
            </div>
          </div>
          <p class="param-note">内存格式示例：512m、4g、128g。这里会写入本次网页运行生成的 config.yaml。</p>
        </details>

        <div class="actions">
          <button onclick="dryRun()">Dry-run</button>
          <button onclick="runWorkflow()">运行流程</button>
          <button class="danger" onclick="stopJob()">停止任务</button>
        </div>
        <p class="muted">网页只选择服务器已有路径，不上传测序数据。</p>
      </section>

      <section class="split">
        <div>
          <div class="section-head">
            <h2>任务进度</h2>
            <span class="badge" id="progressText">0%</span>
          </div>
          <div class="dashboard">
            <div class="metric"><div class="value" id="sampleCount">0</div><div class="label">样本数</div></div>
            <div class="metric"><div class="value" id="chromCount">0</div><div class="label">染色体数</div></div>
            <div class="metric"><div class="value" id="doneJobs">0</div><div class="label">完成任务</div></div>
            <div class="metric"><div class="value" id="totalJobs">0</div><div class="label">总任务</div></div>
          </div>
          <div class="progress-wrap">
            <div class="progress-header">
              <strong id="stageName">等待提交</strong>
              <span class="muted" id="progressDetail">0 / 0 steps</span>
            </div>
            <div class="progress-bar"><div class="progress-fill" id="progressFill"></div></div>
            <div class="stage-list" id="stageList"></div>
          </div>

          <div class="report-actions" id="reportActions">
            <p class="muted">报告已生成，可以直接在网页中查看。</p>
            <button onclick="openReport()">查看报告</button>
          </div>

          <label>当前任务 ID</label>
          <input id="jobId" readonly>
          <input id="configPath" type="hidden">

          <h3>识别样本</h3>
          <div class="table-wrap">
            <table id="sampleTable">
              <thead><tr><th>样本</th><th>FASTQ 前缀</th></tr></thead>
              <tbody></tbody>
            </table>
          </div>
        </div>
        <div>
          <details>
            <summary>查看 Snakemake 原始日志</summary>
            <pre id="log" class="log">等待提交任务。</pre>
          </details>
        </div>
      </section>
    </main>
  </div>

  <div class="modal" id="browserModal">
    <div class="browser">
      <div class="browser-head">
        <h2 id="browserTitle">选择路径</h2>
        <button class="ghost" onclick="closeBrowser()">关闭</button>
      </div>
      <div class="browser-path">
        <input id="browserPath">
        <button class="secondary" onclick="loadBrowser(document.getElementById('browserPath').value)">打开</button>
      </div>
      <div class="browser-list" id="browserList"></div>
      <div class="browser-foot">
        <span class="muted" id="browserHint">选择服务器端文件或目录</span>
        <div class="actions" style="margin:0">
          <button class="secondary" onclick="selectCurrentDirectory()">使用当前目录</button>
          <button onclick="confirmSelection()">确认选择</button>
        </div>
      </div>
    </div>
  </div>

  <script>
    let activeJob = "";
    let timer = null;
    let browserTarget = "";
    let browserMode = "file";
    let browserExt = "";
    let browserSelected = "";
    let browserCurrent = "";
    const stages = [
      ["precheck", "运行前检查"],
      ["qc", "原始质控"],
      ["clean_reads", "数据清洗"],
      ["bwa_mem", "序列比对"],
      ["bam_rmdup", "去重复"],
      ["haplo", "染色体 calling"],
      ["joint_calling", "联合分型"],
      ["snp_filter", "变异过滤"],
      ["vcf_missing", "缺失率统计"],
      ["vcf_to", "格式转换"],
      ["imputation", "基因型填充"],
      ["cnv", "CNV 检测"],
      ["vcf2pca", "PCA"],
      ["vcf2dis", "遗传距离"],
      ["snpeff", "功能注释"],
      ["report", "报告生成"]
    ];

    async function api(path, options = {}) {
      const res = await fetch(path, options);
      const data = await res.json();
      if (!res.ok || data.error) throw new Error(data.error || "request failed");
      return data;
    }

    function lines(id) {
      return document.getElementById(id).value.split(/\r?\n/).map(x => x.trim()).filter(Boolean);
    }

    function numberValue(id, fallback) {
      const value = Number(document.getElementById(id).value || fallback);
      if (!Number.isFinite(value) || value < 1) return fallback;
      return Math.floor(value);
    }

    function textValue(id, fallback) {
      const value = document.getElementById(id).value.trim();
      return value || fallback;
    }

    function updateCounts() {
      document.getElementById("sampleCount").textContent = lines("samples").length;
      document.getElementById("chromCount").textContent = lines("chromosomes").length;
    }

    function collectConfig(dryRunMode) {
      const sampleLines = lines("samples");
      const samples = {};
      for (const line of sampleLines) {
        const parts = line.split(/\t|,/).map(x => x.trim()).filter(Boolean);
        if (parts.length >= 2) samples[parts[0]] = parts[1];
      }
      return {
        fastq_dir: document.getElementById("fastqDir").value.trim(),
        reference: document.getElementById("reference").value.trim(),
        annotation_file: document.getElementById("annotationFile").value.trim(),
        container_image: document.getElementById("containerImage").value.trim(),
        snakemake_cmd: document.getElementById("snakemakeCmd").value.trim(),
        cores: Number(document.getElementById("cores").value || 1),
        singularity_args: document.getElementById("singularityArgs").value.trim(),
        samples: samples,
        chromosomes: lines("chromosomes"),
        modules: {
          vcf2pca: document.getElementById("runPca").checked,
          vcf2dis: document.getElementById("runDis").checked,
          snpeff: document.getElementById("runSnpeff").checked,
          snpeff_indel: document.getElementById("runIndelAnno").checked,
          imputation: document.getElementById("runImputation").checked,
          cnv: document.getElementById("runCnv").checked
        },
        advanced: {
          threads: {
            qc: numberValue("qcThreads", 8),
            fastp: numberValue("fastpThreads", 8),
            bwa: numberValue("bwaThreads", 8),
            pca: numberValue("pcaThreads", 4),
            imputation: numberValue("imputationThreads", 4),
            genomicsdb: numberValue("genomicsdbThreads", 4),
            genomicsdb_batch_size: numberValue("genomicsdbBatchSize", 50),
            cnv_bin_size: numberValue("cnvBinSize", 100),
            cnv_min_len: numberValue("cnvMinLen", 1000)
          },
          memory: {
            haplo_xms: textValue("haploXms", "512m"),
            haplo_xmx: textValue("haploXmx", "4g"),
            genotype_xms: textValue("genotypeXms", "512m"),
            genotype_xmx: textValue("genotypeXmx", "128g"),
            combine_xms: textValue("combineXms", "512m"),
            combine_xmx: textValue("combineXmx", "128g"),
            genomicsdb_xms: textValue("genomicsdbXms", "1g"),
            genomicsdb_xmx: textValue("genomicsdbXmx", "16g"),
            snpeff_xmx: textValue("snpeffXmx", "32g"),
            beagle_xmx: textValue("beagleXmx", "4g")
          },
          cnv: {
            max_evalue: textValue("cnvMaxEvalue", "0.01"),
            max_q0: textValue("cnvMaxQ0", "0.5")
          }
        },
        dry_run: dryRunMode
      };
    }

    function renderSamples(samples) {
      const tbody = document.querySelector("#sampleTable tbody");
      tbody.innerHTML = "";
      for (const [name, prefix] of Object.entries(samples)) {
        const tr = document.createElement("tr");
        tr.innerHTML = `<td>${escapeHtml(name)}</td><td>${escapeHtml(prefix)}</td>`;
        tbody.appendChild(tr);
      }
      updateCounts();
    }

    function escapeHtml(text) {
      return String(text).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
    }

    async function inferSamples() {
      try {
        const fastqDir = encodeURIComponent(document.getElementById("fastqDir").value.trim());
        const data = await api(`/api/infer_samples?dir=${fastqDir}`);
        const text = Object.entries(data.samples).map(([k, v]) => `${k}\t${v}`).join("\n");
        document.getElementById("samples").value = text;
        renderSamples(data.samples);
        document.getElementById("serverStatus").textContent = `samples: ${Object.keys(data.samples).length}`;
      } catch (err) {
        alert(err.message);
      }
    }

    async function inferChromosomes() {
      try {
        const ref = encodeURIComponent(document.getElementById("reference").value.trim());
        const data = await api(`/api/chromosomes?reference=${ref}`);
        document.getElementById("chromosomes").value = data.chromosomes.join("\n");
        updateCounts();
        document.getElementById("serverStatus").textContent = `chromosomes: ${data.chromosomes.length}`;
      } catch (err) {
        alert(err.message);
      }
    }

    async function submit(dryRunMode) {
      try {
        const data = await api("/api/run", {
          method: "POST",
          headers: {"Content-Type": "application/json"},
          body: JSON.stringify(collectConfig(dryRunMode))
        });
        activeJob = data.job_id;
        document.getElementById("jobId").value = data.job_id;
        document.getElementById("configPath").value = data.config_path;
        document.getElementById("jobStatus").textContent = `job: ${data.status}`;
        setProgress({percent: 0, done: 0, total: 0, current_rule: "queued", message: "任务已提交", stages: {}});
        poll();
        if (timer) clearInterval(timer);
        timer = setInterval(poll, 2500);
      } catch (err) {
        alert(err.message);
      }
    }

    function dryRun() { submit(true); }
    function runWorkflow() { submit(false); }

    async function poll() {
      if (!activeJob) return;
      try {
        const data = await api(`/api/status?job_id=${encodeURIComponent(activeJob)}`);
        document.getElementById("jobStatus").textContent = `job: ${data.status}`;
        document.getElementById("log").textContent = data.log || "";
        setProgress(data.progress || {});
        setReportLink(Boolean(data.report_exists) && data.status === "success");
        const logBox = document.getElementById("log");
        logBox.scrollTop = logBox.scrollHeight;
        if (["success", "failed", "stopped"].includes(data.status) && timer) {
          clearInterval(timer);
          timer = null;
        }
      } catch (err) {
        document.getElementById("jobStatus").textContent = "job: status error";
      }
    }

    function setProgress(progress) {
      const percent = Math.max(0, Math.min(100, Math.round(progress.percent || 0)));
      document.getElementById("progressFill").style.width = `${percent}%`;
      document.getElementById("progressText").textContent = `${percent}%`;
      document.getElementById("doneJobs").textContent = progress.done || 0;
      document.getElementById("totalJobs").textContent = progress.total || 0;
      document.getElementById("progressDetail").textContent = `${progress.done || 0} / ${progress.total || 0} steps`;
      document.getElementById("stageName").textContent = progress.message || progress.current_rule || "等待提交";
      renderStages(progress.stages || {}, progress.current_rule || "");
    }

    function setReportLink(visible) {
      document.getElementById("reportActions").classList.toggle("visible", visible);
    }

    function openReport() {
      window.open("/report/ParaChrSNP_report.html", "_blank", "noopener");
    }

    function stageLabel(state) {
      if (state === "completed") return "已完成";
      if (state === "running") return "正在运行";
      return "待运行";
    }

    function renderStages(stageStates, currentRule) {
      const box = document.getElementById("stageList");
      box.innerHTML = "";
      for (const [key, name] of stages) {
        const item = document.createElement("div");
        let state = stageStates[key] || "pending";
        if (currentRule && currentRule.includes(key) && state === "pending") state = "running";
        item.className = `stage ${state}`;
        item.innerHTML = `<div class="name">${escapeHtml(name)}</div><div class="state">${stageLabel(state)}</div>`;
        box.appendChild(item);
      }
    }

    async function openBrowser(target, mode, ext = "") {
      browserTarget = target;
      browserMode = mode;
      browserExt = ext;
      browserSelected = "";
      const currentValue = document.getElementById(target).value.trim();
      document.getElementById("browserTitle").textContent = mode === "dir" ? "选择目录" : "选择文件";
      document.getElementById("browserHint").textContent = mode === "dir" ? "双击目录进入，确认后使用当前目录或选中的目录。" : "双击目录进入，选择文件后确认。";
      document.getElementById("browserModal").classList.add("open");
      await loadBrowser(currentValue || "");
    }

    function closeBrowser() {
      document.getElementById("browserModal").classList.remove("open");
    }

    async function loadBrowser(path) {
      try {
        const url = `/api/browse?path=${encodeURIComponent(path || "")}&mode=${encodeURIComponent(browserMode)}&ext=${encodeURIComponent(browserExt)}`;
        const data = await api(url);
        browserCurrent = data.path || "";
        browserSelected = "";
        document.getElementById("browserPath").value = browserCurrent;
        renderBrowser(data);
      } catch (err) {
        alert(err.message);
      }
    }

    function renderBrowser(data) {
      const list = document.getElementById("browserList");
      list.innerHTML = "";
      for (const item of data.entries) {
        const row = document.createElement("div");
        row.className = "file-row";
        row.dataset.path = item.path;
        row.dataset.kind = item.kind;
        row.innerHTML = `
          <div class="file-kind">${item.kind === "dir" ? "DIR" : "FILE"}</div>
          <div class="file-name">${escapeHtml(item.name)}</div>
          <div class="file-size">${escapeHtml(item.size || "")}</div>
          <div class="file-time">${escapeHtml(item.modified || "")}</div>
        `;
        row.onclick = () => {
          document.querySelectorAll(".file-row").forEach(x => x.classList.remove("selected"));
          row.classList.add("selected");
          browserSelected = item.path;
        };
        row.ondblclick = () => {
          if (item.kind === "dir") loadBrowser(item.path);
          else {
            browserSelected = item.path;
            confirmSelection();
          }
        };
        list.appendChild(row);
      }
    }

    function selectCurrentDirectory() {
      browserSelected = browserCurrent;
      confirmSelection();
    }

    function confirmSelection() {
      const value = browserSelected || browserCurrent;
      if (!value) return;
      document.getElementById(browserTarget).value = value;
      closeBrowser();
      if (browserTarget === "fastqDir") inferSamples();
      if (browserTarget === "reference") inferChromosomes();
    }

    async function stopJob() {
      if (!activeJob) return;
      try {
        await api("/api/stop", {
          method: "POST",
          headers: {"Content-Type": "application/json"},
          body: JSON.stringify({job_id: activeJob})
        });
        poll();
      } catch (err) {
        alert(err.message);
      }
    }

    api("/api/health").then(data => {
      document.getElementById("serverStatus").textContent = `server: ${data.status}`;
    }).catch(() => {
      document.getElementById("serverStatus").textContent = "server: error";
    });
    api("/api/report_status").then(data => {
      setReportLink(Boolean(data.report_exists));
    }).catch(() => {});
    ["samples", "chromosomes"].forEach(id => document.getElementById(id).addEventListener("input", updateCounts));
    renderStages({}, "");
    updateCounts();
  </script>
</body>
</html>
"""


def parse_args():
    """解析命令行参数，并提供英文帮助文档。"""
    parser = argparse.ArgumentParser(
        description="Launch a lightweight web interface for the ParaChrSNP Snakemake workflow.",
        epilog=f"Author: {AUTHOR}",
    )
    parser.add_argument(
        "--host",
        default="0.0.0.0",
        help="Host address used by the web server. Default: 0.0.0.0.",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8088,
        help="Port used by the web server. Default: 8088.",
    )
    parser.add_argument(
        "--allowed-root",
        action="append",
        default=[],
        help="Allowed server-side data root. Can be specified multiple times.",
    )
    return parser.parse_args()


def write_script_log():
    """把脚本信息记录到用户指定的脚本日志中。"""
    log_path = Path("/home/majunpeng/script_log.txt")
    line = (
        f"{time.strftime('%Y-%m-%d %H:%M:%S')}\t{SCRIPT_NAME}\t"
        f"{SCRIPT_FUNCTION}\t{Path(__file__).resolve()}\n"
    )
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(line)


def allowed_roots(extra_roots=None):
    """返回允许网页访问的服务器端目录。"""
    roots = list(DEFAULT_ALLOWED_ROOTS)
    env_roots = os.environ.get("PARACHRSNP_ALLOWED_ROOTS", "")
    for item in env_roots.split(":"):
        if item:
            roots.append(Path(item))
    for item in extra_roots or []:
        roots.append(Path(item))
    valid = []
    for root in roots:
        try:
            if root.exists():
                valid.append(root.resolve())
        except OSError:
            continue
    return valid


def resolve_allowed(path_text, roots):
    """解析用户提供的路径，并限制其必须位于允许访问的目录下。"""
    path = Path(path_text).expanduser()
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    resolved = path.resolve()
    for root in roots:
        if resolved == root or root in resolved.parents:
            return resolved
    raise ValueError(f"Path is outside allowed roots: {path_text}")


def is_allowed_path(path, roots):
    """判断路径是否位于允许访问的目录中。"""
    try:
        resolved = Path(path).resolve()
    except OSError:
        return False
    for root in roots:
        if resolved == root or root in resolved.parents:
            return True
    return False


def format_size(size):
    """把文件大小转换为便于网页显示的字符串。"""
    value = float(size)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if value < 1024 or unit == "TB":
            if unit == "B":
                return f"{int(value)} {unit}"
            return f"{value:.1f} {unit}"
        value /= 1024
    return f"{size} B"


def matches_extension(path, ext_text):
    """根据扩展名过滤文件浏览器中展示的文件。"""
    if not ext_text:
        return True
    name = path.name.lower()
    extensions = [item.strip().lower() for item in ext_text.split(",") if item.strip()]
    return any(name.endswith(ext) for ext in extensions)


def browse_path(path_text, roots, mode="file", ext_text=""):
    """列出允许目录中的文件和子目录，供网页文件浏览器使用。"""
    if not path_text:
        entries = []
        for root in roots:
            entries.append({
                "name": str(root),
                "path": str(root),
                "kind": "dir",
                "size": "",
                "modified": "",
            })
        return {"path": "", "entries": entries}

    path = resolve_allowed(path_text, roots)
    if path.is_file():
        path = path.parent
    if not path.is_dir():
        raise ValueError(f"Directory does not exist: {path}")

    entries = []
    parent = path.parent
    if is_allowed_path(parent, roots) and parent != path:
        entries.append({
            "name": "..",
            "path": str(parent),
            "kind": "dir",
            "size": "",
            "modified": "",
        })

    try:
        children = sorted(
            path.iterdir(),
            key=lambda item: (not item.is_dir(), item.name.lower()),
        )
    except PermissionError:
        raise ValueError(f"No permission to read directory: {path}")

    for child in children:
        try:
            stat = child.stat()
        except OSError:
            continue
        if child.is_dir():
            entries.append({
                "name": child.name,
                "path": str(child.resolve()),
                "kind": "dir",
                "size": "",
                "modified": time.strftime("%Y-%m-%d %H:%M", time.localtime(stat.st_mtime)),
            })
        elif mode != "dir" and matches_extension(child, ext_text):
            entries.append({
                "name": child.name,
                "path": str(child.resolve()),
                "kind": "file",
                "size": format_size(stat.st_size),
                "modified": time.strftime("%Y-%m-%d %H:%M", time.localtime(stat.st_mtime)),
            })
    return {"path": str(path), "entries": entries}


def stage_for_rule(rule_name):
    """根据 rule 名称判断它属于哪个流程阶段。"""
    for stage, rule_names in STAGE_RULES.items():
        if rule_name in rule_names:
            return stage
    return ""


def parse_rule_jobs(log_text):
    """从 Snakemake 日志中解析 jobid 与 rule 的对应关系。"""
    jobs = {}
    pattern = re.compile(
        r"(?m)^(?:localrule|rule)\s+([A-Za-z0-9_]+):(?P<body>.*?)(?=^\[|\Z)",
        re.S,
    )
    for match in pattern.finditer(log_text):
        rule_name = match.group(1)
        job_match = re.search(r"(?m)^\s*jobid:\s*(\d+)\s*$", match.group("body"))
        if job_match:
            jobs[job_match.group(1)] = rule_name
    return jobs


def parse_job_stats_rules(log_text):
    """从 Job stats 表格中解析本次 DAG 计划运行的 rule。"""
    rules = set()
    in_stats = False
    for line in log_text.splitlines():
        if line.startswith("Job stats:"):
            in_stats = True
            continue
        if not in_stats:
            continue
        stripped = line.strip()
        if not stripped:
            if rules:
                break
            continue
        if stripped.startswith("job ") or stripped.startswith("-"):
            continue
        parts = stripped.split()
        if len(parts) >= 2 and parts[0] != "total":
            rules.add(parts[0])
    return rules


def parse_stage_states(log_text, status):
    """根据 Snakemake 日志为每个流程阶段生成 pending/running/completed 状态。"""
    stage_states = {stage: "pending" for stage in STAGE_RULES}
    if not log_text:
        return stage_states

    job_rules = parse_rule_jobs(log_text)
    planned_rules = parse_job_stats_rules(log_text) or set(job_rules.values())
    finished_job_ids = set(re.findall(r"Finished job\s+(\d+)", log_text))
    finished_rules = {job_rules[job_id] for job_id in finished_job_ids if job_id in job_rules}

    for stage, rule_names in STAGE_RULES.items():
        rules = set(rule_names)
        if planned_rules and not rules.intersection(planned_rules):
            stage_states[stage] = "completed"
        if rules.intersection(finished_rules):
            stage_states[stage] = "completed"

    if status == "running":
        for job_id, rule_name in job_rules.items():
            if job_id in finished_job_ids:
                continue
            stage = stage_for_rule(rule_name)
            if stage and stage_states.get(stage) != "completed":
                stage_states[stage] = "running"

    if status == "success" and "This was a dry-run" not in log_text:
        for stage in list(stage_states):
            stage_states[stage] = "completed"
    return stage_states


def parse_progress(log_text, status):
    """从 Snakemake 日志中提取粗略进度信息。"""
    progress = {
        "done": 0,
        "total": 0,
        "percent": 0,
        "current_rule": "",
        "message": "等待提交",
        "stages": parse_stage_states(log_text, status),
    }
    if not log_text:
        if status == "running":
            progress["message"] = "任务启动中"
        return progress

    step_matches = re.findall(r"(\d+)\s+of\s+(\d+)\s+steps\s+\(\d+%\)\s+done", log_text)
    if step_matches:
        done, total = map(int, step_matches[-1])
        progress["done"] = done
        progress["total"] = total
    else:
        total_matches = re.findall(r"(?m)^\s*total\s+(\d+)(?:\s|$)", log_text)
        finished = len(re.findall(r"Finished job \d+", log_text))
        if total_matches:
            progress["total"] = int(total_matches[-1])
        progress["done"] = finished

    if progress["total"] > 0:
        progress["percent"] = round(progress["done"] * 100 / progress["total"])

    rule_matches = re.findall(r"rule\s+([A-Za-z0-9_]+):", log_text)
    if rule_matches:
        progress["current_rule"] = rule_matches[-1]
        progress["message"] = f"当前步骤：{rule_matches[-1]}"

    if "Building DAG of jobs" in log_text and not progress["current_rule"]:
        progress["message"] = "正在构建任务 DAG"
    if "Nothing to be done" in log_text:
        progress.update({"done": 1, "total": 1, "percent": 100, "message": "没有需要运行的任务"})
    if "WorkflowError:" in log_text:
        workflow_error = log_text.split("WorkflowError:", 1)[1].strip().splitlines()
        if workflow_error:
            progress["message"] = f"Snakemake 错误：{workflow_error[0].strip()}"
    if status == "success":
        if progress["total"] > 0:
            progress["done"] = progress["total"]
        progress["percent"] = 100
        progress["message"] = "Dry-run 完成" if "This was a dry-run" in log_text else "任务完成"
    elif status == "failed":
        if progress["message"].startswith("Snakemake 错误"):
            pass
        else:
            progress["message"] = "任务失败"
    elif status == "stopped":
        progress["message"] = "任务已停止"
    elif status == "queued":
        progress["message"] = "任务排队中"
    return progress


def read_json_body(handler):
    """读取 HTTP POST 请求中的 JSON 数据。"""
    length = int(handler.headers.get("Content-Length", "0"))
    raw = handler.rfile.read(length).decode("utf-8")
    return json.loads(raw or "{}")


def send_json(handler, data, status=200):
    """返回 JSON 响应。"""
    body = json.dumps(data, ensure_ascii=False, indent=2).encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", "application/json; charset=utf-8")
    handler.send_header("Content-Length", str(len(body)))
    handler.end_headers()
    handler.wfile.write(body)


def send_text(handler, text, content_type="text/html; charset=utf-8"):
    """返回文本响应。"""
    body = text.encode("utf-8")
    handler.send_response(200)
    handler.send_header("Content-Type", content_type)
    handler.send_header("Content-Length", str(len(body)))
    handler.end_headers()
    handler.wfile.write(body)


def send_file(handler, path, content_type=None):
    """返回项目目录中的只读静态文件。"""
    body = path.read_bytes()
    guessed_type = content_type or mimetypes.guess_type(str(path))[0] or "application/octet-stream"
    handler.send_response(200)
    handler.send_header("Content-Type", guessed_type)
    handler.send_header("Content-Length", str(len(body)))
    handler.end_headers()
    handler.wfile.write(body)


def infer_samples(fastq_dir):
    """从 FASTQ 目录中识别 {sample}.1.fq.gz 和 {sample}.2.fq.gz 样本配对。"""
    samples = {}
    pattern = re.compile(r"(.+)\.1\.(?:f(?:ast)?q)\.gz$")
    for path in sorted(fastq_dir.iterdir()):
        match = pattern.match(path.name)
        if not match:
            continue
        sample = match.group(1)
        mate = fastq_dir / f"{sample}.2.{path.name.rsplit('.1.', 1)[1]}"
        if mate.exists():
            prefix = str((fastq_dir / sample).resolve())
            samples[sample] = prefix
    return samples


def read_chromosomes(reference):
    """读取 FASTA header 中的序列 ID，作为染色体列表。"""
    opener = gzip.open if reference.suffix == ".gz" else open
    chromosomes = []
    with opener(reference, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                chromosomes.append(line[1:].strip().split()[0])
    return chromosomes


def int_param(value, default, minimum=1):
    """把网页提交的线程数转换为正整数，非法值回退到默认值。"""
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return default
    if parsed < minimum:
        return default
    return parsed


def text_param(value, default):
    """读取网页文本参数，空值回退到默认值。"""
    if value is None:
        return default
    value = str(value).strip()
    return value or default


def gatk_java_options(xms, xmx):
    """生成 GATK 命令需要的 --java-options 参数。"""
    return f'--java-options "-Xms{xms} -Xmx{xmx}"'


def xmx_option(value, default):
    """生成普通 Java 程序使用的 -Xmx 参数，兼容用户直接填写 -Xmx4g。"""
    value = text_param(value, default)
    if value.startswith("-Xmx"):
        return value
    return f"-Xmx{value}"


def base_config(payload):
    """根据网页提交参数生成 ParaChrSNP 的 config.yaml 内容。"""
    reference = payload["reference"]
    annotation_file = payload.get("annotation_file", "")
    container_image = payload.get("container_image", "ParaChrSNP.sif")
    modules = payload.get("modules", {})
    advanced = payload.get("advanced", {})
    thread_cfg = advanced.get("threads", {})
    memory_cfg = advanced.get("memory", {})
    cnv_cfg = advanced.get("cnv", {})

    qc_threads = int_param(thread_cfg.get("qc"), 8)
    fastp_threads = int_param(thread_cfg.get("fastp"), 8)
    bwa_threads = int_param(thread_cfg.get("bwa"), 8)
    pca_threads = int_param(thread_cfg.get("pca"), 4)
    imputation_threads = int_param(thread_cfg.get("imputation"), 4)
    genomicsdb_threads = int_param(thread_cfg.get("genomicsdb"), 4)
    genomicsdb_batch_size = int_param(thread_cfg.get("genomicsdb_batch_size"), 50)
    cnv_bin_size = int_param(thread_cfg.get("cnv_bin_size"), 100)
    cnv_min_len = int_param(thread_cfg.get("cnv_min_len"), 1000)

    haplo_java = gatk_java_options(
        text_param(memory_cfg.get("haplo_xms"), "512m"),
        text_param(memory_cfg.get("haplo_xmx"), "4g"),
    )
    genotype_java = gatk_java_options(
        text_param(memory_cfg.get("genotype_xms"), "512m"),
        text_param(memory_cfg.get("genotype_xmx"), "128g"),
    )
    combine_java = gatk_java_options(
        text_param(memory_cfg.get("combine_xms"), "512m"),
        text_param(memory_cfg.get("combine_xmx"), "128g"),
    )
    genomicsdb_java = gatk_java_options(
        text_param(memory_cfg.get("genomicsdb_xms"), "1g"),
        text_param(memory_cfg.get("genomicsdb_xmx"), "16g"),
    )
    snpeff_java = xmx_option(memory_cfg.get("snpeff_xmx"), "32g")
    beagle_java = xmx_option(memory_cfg.get("beagle_xmx"), "4g")

    config = {
        "reference": reference,
        "container": {"image": container_image},
        "samples": payload.get("samples", {}),
        "chromosomes": payload.get("chromosomes", []),
        "params": {
            "qc": {"threads": qc_threads, "extra": "--nozip"},
            "bwa": {"bwa_threads": bwa_threads},
            "fastp": {"fastp_threads": fastp_threads},
            "haplotype_caller": {"java_options": haplo_java},
            "genotype_gvcfs": {"java_options": genotype_java},
            "combine_vcf": {"java_options": combine_java},
            "joint_calling": {
                "method": "genomicsdb",
                "reader_threads": genomicsdb_threads,
                "batch_size": genomicsdb_batch_size,
                "import_java_options": genomicsdb_java,
                "genotype_java_options": genotype_java,
                "gather_java_options": combine_java,
                "extra_import": "",
                "extra_genotype": "",
            },
            "snp_filter": {
                "filter_name": "SNP_filter",
                "thresholds": {
                    "qd": 2.0,
                    "mq": 40.0,
                    "fs": 60.0,
                    "sor": 3.0,
                    "mq_rank_sum": -12.5,
                    "read_pos_rank_sum": -8.0,
                },
            },
            "indel_filter": {
                "filter_name": "INDEL_filter",
                "thresholds": {
                    "qd": 2.0,
                    "fs": 200.0,
                    "sor": 10.0,
                    "mq_rank_sum": -12.5,
                    "read_pos_rank_sum": -8.0,
                },
            },
            "vcf_missing": {
                "output_prefix": "missing/combined.snp.filtered",
                "plot_prefix": "missing/combined.snp.filtered.miss_check",
                "plot_script": "scripts/plink-missing.py",
                "extra": "--allow-extra-chr",
            },
            "vcf_convert": {
                "input_vcf": "result_vcfs/combined.snp.filtered.vcf.gz",
                "output_prefix": "format_convert/combined.snp.filtered",
                "plink_extra": "--allow-extra-chr",
                "tassel_executable": "run_pipeline.pl",
                "tassel_memory": "-Xmx10g",
                "tassel_extra": "",
            },
            "imputation": {
                "enabled": bool(modules.get("imputation", False)),
                "input_vcf": "result_vcfs/combined.snp.filtered.vcf.gz",
                "output_prefix": "imputation/combined.snp.filtered.beagle",
                "jar": "",
                "java_options": beagle_java,
                "threads": imputation_threads,
                "extra": "",
            },
            "cnv": {
                "enabled": bool(modules.get("cnv", False)),
                "software": "cnvnator",
                "executable": "cnvnator",
                "bin_size": cnv_bin_size,
                "reference_dir": str(Path(reference).parent) if str(Path(reference).parent) != "." else ".",
                "min_len": cnv_min_len,
                "max_evalue": text_param(cnv_cfg.get("max_evalue"), "0.01"),
                "max_q0": text_param(cnv_cfg.get("max_q0"), "0.5"),
                "extra": "",
            },
            "vcf2pca": {
                "enabled": bool(modules.get("vcf2pca", True)),
                "executable": "VCF2PCACluster",
                "sample_group": "",
                "output_prefix": "pca/ParaChrSNP",
                "threads": pca_threads,
                "extra": "",
            },
            "vcf2dis": {
                "enabled": bool(modules.get("vcf2dis", True)),
                "executable": "VCF2Dis",
                "sample_group": "",
                "output_matrix": "dis/ParaChrSNP.p_dis.mat",
                "output_tree": "dis/ParaChrSNP.p_dis.nwk",
                "tree_method": 1,
                "extra": "",
            },
            "snpeff": {
                "enabled": bool(modules.get("snpeff", False)),
                "annotate_snp": True,
                "annotate_indel": bool(modules.get("snpeff_indel", False)),
                "executable": "snpEff",
                "genome_name": "custom_genome",
                "data_dir": "annotation/snpeff_data",
                "config_file": "annotation/snpeff.config",
                "genome_fasta": reference,
                "annotation_file": annotation_file,
                "annotation_format": "gff3",
                "output_prefix": "annotation/combined",
                "database_done": "annotation/snpeff_db.done",
                "java_options": snpeff_java,
                "build_check_options": "-noCheckCds -noCheckProtein",
                "extra": "",
            },
        },
    }
    return config


def write_config(payload, job_dir):
    """写出本次任务使用的 config.yaml。"""
    if yaml is None:
        raise RuntimeError("PyYAML is required to write config.yaml.")
    config = base_config(payload)
    config_path = job_dir / "config.yaml"
    with config_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(config, handle, allow_unicode=True, sort_keys=False)
    return config_path


def run_process(job_id, command, log_path):
    """在后台运行 Snakemake，并实时记录任务状态。"""
    env = os.environ.copy()
    wrapper_dir = PROJECT_ROOT / "web" / "bin"
    if wrapper_dir.is_dir():
        env["PATH"] = f"{wrapper_dir}{os.pathsep}{env.get('PATH', '')}"
    with JOBS_LOCK:
        JOBS[job_id]["status"] = "running"
        JOBS[job_id]["command"] = command
    with log_path.open("a", encoding="utf-8") as log:
        log.write("$ " + " ".join(command) + "\n\n")
        log.flush()
        try:
            process = subprocess.Popen(
                command,
                cwd=str(PROJECT_ROOT),
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                env=env,
                preexec_fn=os.setsid,
            )
        except FileNotFoundError as exc:
            log.write(f"\nERROR: command not found: {command[0]}\n")
            log.write("Please activate the Snakemake environment before starting this web server, ")
            log.write("or set the Snakemake command field to an absolute executable path.\n")
            log.flush()
            with JOBS_LOCK:
                JOBS[job_id]["status"] = "failed"
                JOBS[job_id]["return_code"] = 127
            return
        with JOBS_LOCK:
            JOBS[job_id]["process"] = process
        return_code = process.wait()
    with JOBS_LOCK:
        if JOBS[job_id].get("status") == "stopped":
            return
        JOBS[job_id]["status"] = "success" if return_code == 0 else "failed"
        JOBS[job_id]["return_code"] = return_code


class ParaChrSNPHandler(BaseHTTPRequestHandler):
    """处理 ParaChrSNP Web Launcher 的 HTTP 请求。"""

    server_version = "ParaChrSNPWeb/0.1"

    def do_GET(self):
        """处理网页和 API GET 请求。"""
        try:
            parsed = urlparse(self.path)
            if parsed.path == "/":
                page = INDEX_HTML.replace("__DEFAULT_SNAKEMAKE__", html.escape(DEFAULT_SNAKEMAKE, quote=True))
                send_text(self, page)
            elif parsed.path == "/icon.png":
                icon = PROJECT_ROOT / "figures" / "ParaChrSNP_icon.png"
                body = icon.read_bytes()
                self.send_response(200)
                self.send_header("Content-Type", "image/png")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
            elif parsed.path == "/report/ParaChrSNP_report.html":
                report = PROJECT_ROOT / "reports" / "ParaChrSNP_report.html"
                if not report.is_file():
                    self.send_error(404, "Report has not been generated.")
                    return
                send_text(self, report.read_text(encoding="utf-8", errors="replace"))
            elif parsed.path == "/report/ParaChrSNP_summary.tsv":
                summary = PROJECT_ROOT / "reports" / "ParaChrSNP_summary.tsv"
                if not summary.is_file():
                    self.send_error(404, "Summary has not been generated.")
                    return
                send_text(self, summary.read_text(encoding="utf-8", errors="replace"), "text/tab-separated-values; charset=utf-8")
            elif parsed.path.startswith("/outputs/"):
                relative = parsed.path.removeprefix("/outputs/")
                target = (PROJECT_ROOT / relative).resolve()
                if PROJECT_ROOT not in target.parents or not target.is_file():
                    self.send_error(404, "Output file is not available.")
                    return
                allowed_suffixes = {
                    ".html", ".htm", ".pdf", ".png", ".svg", ".txt", ".tsv",
                    ".csv", ".nwk", ".mat", ".log",
                }
                if target.suffix.lower() not in allowed_suffixes:
                    self.send_error(403, "This output file type is not allowed.")
                    return
                send_file(self, target)
            elif parsed.path == "/api/health":
                send_json(self, {"status": "ready", "project_root": str(PROJECT_ROOT)})
            elif parsed.path == "/api/report_status":
                report = PROJECT_ROOT / "reports" / "ParaChrSNP_report.html"
                send_json(
                    self,
                    {
                        "report_exists": report.is_file(),
                        "report_url": "/report/ParaChrSNP_report.html",
                    },
                )
            elif parsed.path == "/api/browse":
                query = parse_qs(parsed.query)
                path = query.get("path", [""])[0]
                mode = query.get("mode", ["file"])[0]
                ext_text = query.get("ext", [""])[0]
                send_json(self, browse_path(path, self.server.allowed_roots, mode, ext_text))
            elif parsed.path == "/api/infer_samples":
                query = parse_qs(parsed.query)
                fastq_dir = resolve_allowed(query.get("dir", [""])[0], self.server.allowed_roots)
                if not fastq_dir.is_dir():
                    raise ValueError(f"FASTQ directory does not exist: {fastq_dir}")
                send_json(self, {"samples": infer_samples(fastq_dir)})
            elif parsed.path == "/api/chromosomes":
                query = parse_qs(parsed.query)
                reference = resolve_allowed(query.get("reference", [""])[0], self.server.allowed_roots)
                if not reference.is_file():
                    raise ValueError(f"Reference FASTA does not exist: {reference}")
                send_json(self, {"chromosomes": read_chromosomes(reference)})
            elif parsed.path == "/api/status":
                query = parse_qs(parsed.query)
                job_id = query.get("job_id", [""])[0]
                with JOBS_LOCK:
                    job = JOBS.get(job_id)
                if not job:
                    raise ValueError("Unknown job_id.")
                log = ""
                if job["log_path"].exists():
                    text = job["log_path"].read_text(encoding="utf-8", errors="replace")
                    log = text[-60000:]
                progress = parse_progress(log, job["status"])
                send_json(
                    self,
                    {
                        "job_id": job_id,
                        "status": job["status"],
                        "config_path": str(job["config_path"]),
                        "progress": progress,
                        "report_exists": (PROJECT_ROOT / "reports" / "ParaChrSNP_report.html").is_file(),
                        "report_url": "/report/ParaChrSNP_report.html",
                        "log": log,
                    },
                )
            else:
                self.send_error(404)
        except Exception as exc:
            send_json(self, {"error": str(exc)}, status=400)

    def do_POST(self):
        """处理任务提交和停止请求。"""
        try:
            parsed = urlparse(self.path)
            if parsed.path == "/api/run":
                payload = read_json_body(self)
                self.validate_payload(payload)
                job_id = time.strftime("%Y%m%d_%H%M%S_") + uuid.uuid4().hex[:8]
                job_dir = RUN_ROOT / job_id
                job_dir.mkdir(parents=True, exist_ok=True)
                config_path = write_config(payload, job_dir)
                log_path = job_dir / "snakemake.log"
                snakemake_cmd = payload.get("snakemake_cmd", "").strip()
                if not snakemake_cmd:
                    snakemake_cmd = os.environ.get("PARACHRSNP_SNAKEMAKE", DEFAULT_SNAKEMAKE)
                command = shlex.split(snakemake_cmd) + [
                    "--snakefile",
                    "Snakefile",
                    "--configfile",
                    str(config_path),
                    "--cores",
                    str(int(payload.get("cores", 1))),
                    "--use-singularity",
                ]
                singularity_args = payload.get("singularity_args", "").strip()
                if singularity_args:
                    command.extend(["--singularity-args", singularity_args])
                if payload.get("dry_run", False):
                    command.append("-n")
                with JOBS_LOCK:
                    JOBS[job_id] = {
                        "status": "queued",
                        "config_path": config_path,
                        "log_path": log_path,
                        "process": None,
                    }
                thread = threading.Thread(target=run_process, args=(job_id, command, log_path), daemon=True)
                thread.start()
                send_json(
                    self,
                    {
                        "job_id": job_id,
                        "status": "queued",
                        "config_path": str(config_path),
                    },
                )
            elif parsed.path == "/api/stop":
                payload = read_json_body(self)
                job_id = payload.get("job_id", "")
                with JOBS_LOCK:
                    job = JOBS.get(job_id)
                    process = job.get("process") if job else None
                    if job:
                        job["status"] = "stopped"
                if not job:
                    raise ValueError("Unknown job_id.")
                if process and process.poll() is None:
                    os.killpg(os.getpgid(process.pid), signal.SIGTERM)
                send_json(self, {"job_id": job_id, "status": "stopped"})
            else:
                self.send_error(404)
        except Exception as exc:
            send_json(self, {"error": str(exc)}, status=400)

    def validate_payload(self, payload):
        """检查网页提交的关键路径和参数是否可用。"""
        resolve_allowed(payload.get("reference", ""), self.server.allowed_roots)
        fastq_dir = resolve_allowed(payload.get("fastq_dir", ""), self.server.allowed_roots)
        if not fastq_dir.is_dir():
            raise ValueError(f"FASTQ directory does not exist: {fastq_dir}")
        annotation_file = payload.get("annotation_file", "").strip()
        if annotation_file:
            resolve_allowed(annotation_file, self.server.allowed_roots)
        samples = payload.get("samples", {})
        if not samples:
            raise ValueError("No samples were provided.")
        chromosomes = payload.get("chromosomes", [])
        if not chromosomes:
            raise ValueError("No chromosomes were provided.")
        cores = int(payload.get("cores", 1))
        if cores < 1:
            raise ValueError("cores must be >= 1.")

    def log_message(self, format_text, *args):
        """减少终端 HTTP 访问日志噪音。"""
        safe_message = html.escape(format_text % args)
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {self.address_string()} {safe_message}")


def main():
    """启动 ParaChrSNP 网页服务。"""
    args = parse_args()
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    try:
        write_script_log()
    except OSError:
        pass
    roots = allowed_roots(args.allowed_root)
    server = ThreadingHTTPServer((args.host, args.port), ParaChrSNPHandler)
    server.allowed_roots = roots
    print(f"ParaChrSNP web launcher: http://{args.host}:{args.port}")
    print("Allowed roots:")
    for root in roots:
        print(f"  - {root}")
    server.serve_forever()


if __name__ == "__main__":
    main()
