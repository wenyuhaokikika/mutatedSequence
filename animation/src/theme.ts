export const COLORS = {
  cream: '#FFF8F0',
  paper: '#FFFDFC',
  ink: '#2D2A32',
  muted: '#756C7A',
  coral: '#FF6F61',
  coralSoft: '#FFE1DB',
  peach: '#FFB7A5',
  teal: '#2A9D8F',
  tealSoft: '#DDF3EE',
  mint: '#B7E4C7',
  lavender: '#CDB4DB',
  lilac: '#F2EAF8',
  gold: '#F4C95D',
  goldSoft: '#FFF2C7',
  green: '#2E7D5B',
  red: '#C44536',
  line: 'rgba(45, 42, 50, 0.12)',
  shadow: '0 24px 60px rgba(91, 55, 68, 0.12)',
  shadowSoft: '0 12px 32px rgba(91, 55, 68, 0.08)',
};

export const FONT =
  '"Noto Sans SC", "Noto Sans CJK SC", "Source Han Sans SC", "Microsoft YaHei", "PingFang SC", sans-serif';
export const MONO =
  '"JetBrains Mono", "Cascadia Code", "SFMono-Regular", Consolas, monospace';

export const card = {
  background: 'rgba(255, 253, 252, 0.88)',
  border: `1px solid ${COLORS.line}`,
  borderRadius: 26,
  boxShadow: COLORS.shadowSoft,
} as const;

export const sequenceGradient = `linear-gradient(90deg, ${COLORS.coral}, ${COLORS.gold})`;
