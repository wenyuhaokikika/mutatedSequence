import React from 'react';
import {AbsoluteFill} from 'remotion';
import {COLORS, FONT, MONO, card} from '../theme';
import {Pill, SceneShell, fadeUp, useEnter} from '../components/ui';
import type {SceneProps} from './types';

const Param: React.FC<{
  flag: string;
  meaning: string;
  delay: number;
}> = ({flag, meaning, delay}) => {
  const p = useEnter(delay);
  return (
    <div
      style={{
        display: 'grid',
        gridTemplateColumns: '190px 1fr',
        gap: 12,
        alignItems: 'baseline',
        padding: '10px 0',
        borderBottom: `1px dashed ${COLORS.line}`,
        opacity: p,
        transform: `translateX(${(1 - p) * 14}px)`,
      }}
    >
      <span style={{fontFamily: MONO, fontSize: 12, fontWeight: 800, color: COLORS.coral}}>{flag}</span>
      <span style={{fontSize: 13, color: COLORS.muted, lineHeight: 1.4}}>{meaning}</span>
    </div>
  );
};

export const DockerScene: React.FC<SceneProps> = ({scene, index, total}) => {
  const title = useEnter(2, 18);
  const cmd = useEnter(12);
  const params = useEnter(26);
  return (
    <SceneShell scene={scene} index={index} total={total} accent={COLORS.teal}>
      <AbsoluteFill>
        <div style={{display: 'flex', alignItems: 'flex-end', justifyContent: 'space-between'}}>
          <div style={fadeUp(title, 18)}>
            <Pill color={COLORS.teal} bg={COLORS.tealSoft} style={{marginBottom: 14}}>
              no setup · docker
            </Pill>
            <h2
              style={{
                margin: 0,
                fontFamily: FONT,
                fontSize: 52,
                lineHeight: 1.15,
                letterSpacing: -2,
                fontWeight: 900,
              }}
            >
              不想装环境？
              <span style={{color: COLORS.teal}}>一行 Docker</span>
            </h2>
          </div>
          <div style={{fontSize: 14, color: COLORS.muted, lineHeight: 1.6, textAlign: 'right'}}>
            镜像只装环境；代码和数据
            <br />
            从挂载目录进来。
          </div>
        </div>

        <div style={{display: 'flex', gap: 28, marginTop: 30, alignItems: 'flex-start'}}>
          <div
            style={{
              width: 960,
              ...card,
              padding: '20px 22px',
              background: '#2D2A32',
              border: '1px solid rgba(255,255,255,0.12)',
              ...fadeUp(cmd, 18),
            }}
          >
            <div style={{display: 'flex', gap: 7, marginBottom: 14}}>
              {[COLORS.coral, COLORS.gold, COLORS.teal].map((c) => (
                <span key={c} style={{width: 9, height: 9, borderRadius: 9, background: c, opacity: 0.8}} />
              ))}
              <span style={{marginLeft: 'auto', fontFamily: MONO, fontSize: 10, color: 'rgba(255,255,255,0.4)'}}>
                terminal
              </span>
            </div>
            <pre
              style={{
                margin: 0,
                fontFamily: MONO,
                fontSize: 15,
                lineHeight: 1.75,
                color: '#F8F0E8',
                whiteSpace: 'pre-wrap',
              }}
            >
              <span style={{color: COLORS.gold}}>docker pull</span>{' '}
              <span style={{color: COLORS.teal}}>wenyuhaobuemoji/mutatedsequence:2.0.0</span>
              {'\n\n'}
              <span style={{color: COLORS.gold}}>docker run</span> -it --rm \{'\n'}
              {'  '}-v <span style={{color: COLORS.coral}}>"$PWD":/workspace</span> \{'\n'}
              {'  '}-w /workspace \{'\n'}
              {'  '}wenyuhaobuemoji/mutatedsequence:2.0.0 \{'\n'}
              {'  '}python -m mutatedseq run \{'\n'}
              {'    '}--genome genome.fa \{'\n'}
              {'    '}--gtf annotation.gtf \{'\n'}
              {'    '}--vcf sample.vcf \{'\n'}
              {'    '}--out out --sample S1
              {'\n'}
              <span style={{color: 'rgba(255,255,255,0.4)'}}>
                # optional: --cds cds.fa.gz --cdna cdna.fa.gz --pep pep.fa.gz
              </span>
            </pre>
          </div>

          <div style={{flex: 1, ...fadeUp(params, 16)}}>
            <div
              style={{
                ...card,
                padding: '14px 18px 8px',
                background: 'rgba(255,253,252,0.9)',
              }}
            >
              <Param flag="docker pull" meaning="把装好环境的镜像拉下来" delay={30} />
              <Param flag="-it --rm" meaning="开一个交互终端，用完自动删容器" delay={34} />
              <Param flag='-v "$PWD":/workspace' meaning="把当前目录挂进容器；代码和数据都在这里" delay={38} />
              <Param flag="-w /workspace" meaning="容器里的工作目录" delay={42} />
              <Param flag="--genome / --gtf / --vcf" meaning="参考基因组、注释、变异表" delay={46} />
              <Param flag="--out" meaning="输出目录" delay={50} />
              <Param flag="--sample" meaning="选择 VCF 里的样本" delay={54} />
              <Param flag="--cds / --cdna / --pep" meaning="可选的校验 FASTA" delay={58} />
            </div>
          </div>
        </div>
      </AbsoluteFill>
    </SceneShell>
  );
};
