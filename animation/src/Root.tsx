import React from 'react';
import {Composition} from 'remotion';
import {Explainer} from './Explainer';
import timing from './timing.json';

export const RemotionRoot: React.FC = () => {
  return (
    <Composition
      id="Explainer"
      component={Explainer}
      durationInFrames={timing.totalFrames}
      fps={timing.fps}
      width={timing.width}
      height={timing.height}
    />
  );
};
