import React from 'react';
import {AbsoluteFill, Audio, Sequence, staticFile} from 'remotion';
import {COLORS} from './theme';
import {SceneRenderer} from './scenes';
import timing from './timing.json';

export const Explainer: React.FC = () => {
  return (
    <AbsoluteFill style={{backgroundColor: COLORS.cream}}>
      {timing.scenes.map((scene, index) => (
        <Sequence
          key={scene.id}
          from={scene.from}
          durationInFrames={scene.durationInFrames}
          name={scene.id}
        >
          <SceneRenderer
            scene={scene}
            index={index}
            total={timing.scenes.length}
          />
          <Audio src={staticFile(scene.audio)} volume={0.96} />
        </Sequence>
      ))}
    </AbsoluteFill>
  );
};
