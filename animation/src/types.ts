export type SceneTiming = {
  id: string;
  title: string;
  narration: string;
  captions: string[];
  audio: string;
  seconds: number;
  durationInFrames: number;
  from: number;
};

export type Timing = {
  fps: number;
  width: number;
  height: number;
  voice: string;
  rate: string;
  pitch: string;
  totalFrames: number;
  totalSeconds: number;
  scenes: SceneTiming[];
};
